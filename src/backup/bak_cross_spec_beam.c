#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "geo_seis.h"

#define PI 3.1415926535

typedef struct para {
    int   nsac;
    float f;
    float s1;
    float s2;
    float b1;
    float b2;
    int   ns;
    int   nb;
    fftw_complex *fft;
    float **pos;
} ARRAYPARA;

void array_spectra(ARRAYPARA p, fftw_complex **sp, char *system);
int  near_pow2(int n);
void no_spa(char *ps);

int main(int argc, char *argv[]) {
    if(argc != 2) {
        fprintf(stderr, "Usage: cross_spec_beam input.para\n");
        exit(1);
    }

    char ss[256], system[128], infile[128], outfile[128];
    int i, j, nfft, fn, index, nsac, n1, n2;
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
                sscanf(ss, "%f", &p.f); break;
            case 5:
                sscanf(ss, "%f", &p.s1); break;
            case 6:
                sscanf(ss, "%f", &p.s2); break;
            case 7:
                sscanf(ss, "%f", &p.b1); break;
            case 8:
                sscanf(ss, "%f", &p.b2); break;
            case 9:
                sscanf(ss, "%d", &p.ns); break;
            case 10:
                sscanf(ss, "%d", &p.nb); break;
            case 11:
                sscanf(ss, "%s", system); break;
            case 12:
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
    n1 = (int) ( (t1-hd.b) / hd.delta );
    n2 = (int) ( (t2-hd.b) / hd.delta );
    nfft  = near_pow2(n2-n1+1);
    fn  = (int) (p.f * nfft * hd.delta);
    fftw_complex *in, *out;
    fftw_plan plan;
    in  = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nfft);
    out = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nfft);

    p.fft = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nsac );
    p.pos  = (float **) malloc( sizeof(float *) * nsac );
    for (i = 0; i < nsac; i ++)
        p.pos[i] = (float *) malloc( sizeof(float) * 2 );

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
        p.fft[index][0] = out[fn][0];
        p.fft[index][1] = out[fn][1];
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
            fprintf(fout, "%g ", sp[i][j][0]*sp[i][j][0]+sp[i][j][1]*sp[i][j][1]);
        if(j == p.nb)
            fprintf(fout, "\n");
    }
    fclose(fout);

// Release dynamic memories.
    printf("Releasing dynamic memories in the main program ...\n");
    for(i = 0; i < nsac; i ++)
        free(p.pos[i]);
    free(p.pos); fftw_free(p.fft);
    for(i = 0; i < p.ns; i ++)
        fftw_free(sp[i]);
    fftw_free(sp);
}

void array_spectra( ARRAYPARA p, fftw_complex **sp, char *system) {
    int i, j, k, nccf;
    float tmp1, tmp2, slow, baz, ds, db, shift;
    float **rpos;
    fftw_complex *ccfd;

    nccf = p.nsac * (p.nsac-1) / 2;
    ccfd = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nccf );
    rpos = (float **) malloc( sizeof(float *) * nccf );
    for(i = 0; i < nccf; i ++)
        rpos[i] = (float *) malloc( sizeof(float) * 2 );

// Initialization of array processing spectra.
    printf("Initialization of array processing spectra ...\n");
    for(i = 0; i < p.ns; i ++)
        for(j = 0; j < p.nb; j ++){
            sp[i][j][0] = 0.;
            sp[i][j][1] = 0.;
    }

// Compute cross-spectra.
    printf("Obtaining cross-spectra ...\n");
    k = 0;
    for(i = 0; i < p.nsac-1; i ++)
        for(j = i+1; j < p.nsac; j ++){
            ccfd[k][0] = p.fft[i][0] * p.fft[j][0] + p.fft[i][1] * p.fft[j][1]; 
            ccfd[k][1] = p.fft[j][0] * p.fft[i][1] - p.fft[i][0] * p.fft[j][1];
            k ++;
    }
// Obtain relative coordinates of seismic array.
    printf("Obtaining relative positions ...\n");
    if( !(strcmp("lonlat", system))) {
        static float g = 111.19492664455;
        float dist, az, res[3];
        k = 0;
        for(i = 0; i < p.nsac-1; i ++)
            for(j = i+1; j < p.nsac; j ++){
                cal_dist_az_baz(p.pos[i][1], p.pos[i][0], \
                                p.pos[j][1], p.pos[j][0], \
                                res);
                dist = res[0] * g;
                az = res[1] / 180 * PI;
                rpos[k][0] = dist * sin(az);
                rpos[k][1] = dist * cos(az);
                k ++;
            }
    } else if( !(strcmp("km", system))) {
        k = 0;
        for(i = 0; i < p.nsac-1; i ++)
            for(j = i+1; j < p.nsac; j ++){
                rpos[k][0] = p.pos[j][0] - p.pos[i][0];
                rpos[k][1] = p.pos[j][1] - p.pos[i][1];
                k ++;
            }
    } else {
        k = 0;
        for(i = 0; i < p.nsac-1; i ++)
            for(j = i+1; j < p.nsac; j ++) {
                rpos[k][0] = (p.pos[j][0] - p.pos[i][0]) / 1e3;
                rpos[k][1] = (p.pos[j][1] - p.pos[i][1]) / 1e3;
                k ++;
            }
    }

    tmp1 = -2.* PI * p.f;
    ds   = (p.s2 - p.s1) / (p.ns - 1);
    db   = (p.b2 - p.b1) / (p.nb - 1);
    printf("Scanning slowness and back-azimuth ...\n");
    for(i = 0; i < p.ns; i ++) {
        slow = p.s1 + i * ds;
        for(j = 0; j < p.nb; j ++){
            baz = (p.b1+j*db-180.) / 180. * PI;
            for(k = 0; k < nccf; k ++){
                shift = slow * (sin(baz)*rpos[k][0]+cos(baz)*rpos[k][1]);
                tmp2 = tmp1 * shift;
                sp[i][j][0] += ( ccfd[k][0] * cos(tmp2) - ccfd[k][1] * sin(tmp2) );
                sp[i][j][1] += ( ccfd[k][1] * cos(tmp2) + ccfd[k][0] * sin(tmp2) );
            }        
        }
    }
    printf("Releasing dynamic memories in array processing ...\n");
// Release dynamic memories.
    for(i = 0; i < nccf; i ++) {
        free(rpos[i]);
    }
    free(rpos); fftw_free(ccfd);
}

int near_pow2(int n) {
    int m, i;
    i = (int) (log((float)n) / log(2.) + 1.);
    m = (int) pow(2., i);
    return m;
}

// De-space.
void no_spa(char *ps) {
    char *pt = ps;
    while ( *ps != '\0' ) {
        if ( *ps != ' ' && *ps != '\n' ) {
            *pt++ = *ps;
        }
        ++ps;
    }
    *pt = '\0';
}
