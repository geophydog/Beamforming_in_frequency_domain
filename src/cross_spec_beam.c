#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "geo_seis.h"

int main(int argc, char *argv[]) {
    if(argc != 2) {
        fprintf(stderr, "Usage: cross_spec_beam input.para\n");
        exit(1);
    }

    char ss[256], system[128], infile[128], outfile[128];
    int i, j, nfft, fn1, fn2, index, nsac, n1, n2;
    float t1, t2, *data;
    fftw_complex **sp;
    FILE *fp, *fin;
    SACHEAD hd;
    ARRAYPARA p;

// Read in parameters.
    fp = fopen(argv[1], "r");
    index = 1;
    while(fgets(ss, 256, fp)) {
        switch(index) {
            case 1:
                sscanf(ss, "%s", infile); break;
            case 2:
                sscanf(ss, "%f", &t1); break;
            case 3:
                sscanf(ss, "%f", &t2); break;
            case 4:
                sscanf(ss, "%f", &p.f1); break;
            case 5:
                sscanf(ss, "%f", &p.f2); break;
            case 6:
                sscanf(ss, "%f", &p.s1); break;
            case 7:
                sscanf(ss, "%f", &p.s2); break;
            case 8:
                sscanf(ss, "%f", &p.b1); break;
            case 9:
                sscanf(ss, "%f", &p.b2); break;
            case 10:
                sscanf(ss, "%d", &p.ns); break;
            case 11:
                sscanf(ss, "%d", &p.nb); break;
            case 12:
                sscanf(ss, "%s", system); break;
            case 13:
                sscanf(ss, "%s", outfile); break;
            default:
                fprintf(stderr, "Ignore the other lines!\n");
        }
        index ++;
    }
    fclose(fp);

// Get the number of SAC files.
    nsac = 0;
    fin = fopen(infile, "r");
    while(fgets(ss, 256, fin)) {
        no_spa(ss);
        nsac ++;
    }
    fclose(fin);

// Allocate dynamic memories of FFT and coordinate arrays.
    data = read_sac(ss, &hd);
    n1   = (int) ( (t1-hd.b) / hd.delta );
    n2   = (int) ( (t2-hd.b) / hd.delta );
    nfft = near_pow2(n2-n1+1);
    fn1  = (int) (p.f1 * nfft * hd.delta);
    fn2  = (int) (p.f2 * nfft * hd.delta);
    p.fn = fn2 - fn1 + 1;

// Check frequency range and Nyquist frequency.
    if((p.fn < 2) || (p.f2 >= 0.5/hd.delta)) {
        fprintf(stderr, "Give larger freuqnecy range or check Nyquist frequency!!!\n");
        exit(1);
    }

    fftw_complex *in, *out;
    fftw_plan plan;
    in  = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nfft);
    out = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nfft);
    
    p.fft = (fftw_complex **) fftw_malloc( sizeof(fftw_complex *) * p.fn );
    p.pos  = (float **) malloc( sizeof(float *) * nsac );
    for (i = 0; i < nsac; i ++)
        p.pos[i] = (float *) malloc( sizeof(float) * 2 );
    for(i = 0; i < p.fn; i ++)
        p.fft[i] = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nsac );

// Obtain Fourier spectra and coordinates of seismic stations.
    printf("Computing Fourier spectra ...\n");    
    index = 0;
    fin = fopen(infile, "r");
    while ( fgets(ss, 256, fin) ) {
        no_spa(ss);
        data = read_sac(ss, &hd);
        n1 = (int) ((t1-hd.b)/hd.delta);
        n2 = (int) ((t2-hd.b)/hd.delta);
        for(i = 0; i < nfft; i ++){
            if(i < (n2-n1+1))
                in[i][0] = data[i+n1];
            else
                in[i][0] = 0.;
            in[i][1] = 0.;
        }
        plan = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        for(i = 0; i < p.fn; i ++) {
            p.fft[i][index][0] = out[fn1+i][0];
            p.fft[i][index][1] = out[fn1+i][1];
        }
        p.pos[index][0] = hd.stlo;
        p.pos[index][1] = hd.stla;
        index ++;
    }
    fftw_destroy_plan(plan);
    fclose(fin);
    fftw_free(in); fftw_free(out);

// Setting parameters of cross-spectral beamforming.
    printf("Setting parameters ...\n");
    p.nsac = nsac;
    sp = (fftw_complex **) fftw_malloc( sizeof(fftw_complex *) * p.ns );
    for(i = 0; i < p.ns; i ++)
        sp[i] = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * p.nb );

// Do cross-spectral beamforming.
    printf("Scanning beamforming spectra ...\n");
    array_spectra(p, sp, system);

// Write out beamforming spectra.
    printf("Writing out beamforming spectra ...\n");
    FILE *fout;
    fout = fopen(outfile, "w");
    for(i = 0; i < p.ns; i ++) {
        for(j = 0; j < p.nb; j ++)
            fprintf(fout, "%g+%gj ", sp[i][j][0], sp[i][j][1]);
        if(j == p.nb)
            fprintf(fout, "\n");
    }
    fclose(fout);

// Release dynamic memories.
    printf("Releasing dynamic memories in the main program ...\n");
    for(i = 0; i < nsac; i ++)
        free(p.pos[i]);
    for(i = 0; i < p.fn; i ++)
        fftw_free(p.fft[i]);
    free(p.pos); fftw_free(p.fft);
    for(i = 0; i < p.ns; i ++)
        fftw_free(sp[i]);
    fftw_free(sp);

    return 0;
}
