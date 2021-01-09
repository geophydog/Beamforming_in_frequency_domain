/*******************************************************************************
 *                                  geo_seis.c                                 *
 *  SAC I/O functions:                                                         *
 *      read_sac         read SAC binary data                                  *
 *      read_head_in     read in SAC header                                    *
 *                       by Dongdong Tian @USTC                                *
 *  GEODETIC functions:                                                        *
 *      cal_dist_az      calculate distance in degree ans azimuth in degree    *
 *                        by Xuping Feng @SUSTech                              *
 *  ARRAY PROCESSING functions:                                                *
 *      cross_spectra_array Computing beam-forming spectra with cross-spectra  *
 *                                                                             *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fftw3.h>
#include "geo_seis.h"

/* a SAC structure containing all null values */
static SACHEAD sac_null = {
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  -12345 , -12345 , -12345 , -12345 , -12345 ,
  { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }
};
/* -------------- function prototype for local use ------------ */
static void    byte_swap       (char *pt, size_t n);
static int     check_sac_nvhdr (const int nvhdr);
static void    map_chdr_in     (char *memar, char *buff);
static int read_head_in(const char *name, SACHEAD *hd, FILE *strm);

/* --------------------------------------------------------------
 *  byte_swap : reverse the byte order of 4 bytes int/float.
 *
 *  IN:
 *      char    *pt : pointer to byte array
 *      size_t   n  : number of bytes
 *  Return: none
 *
 *  Notes:
 *      For 4 bytes,
 *      byte swapping means taking [0][1][2][3],
 *      and turning it into [3][2][1][0]
  -------------------------------------------------------------- */
static void byte_swap(char *pt, size_t n)
{
    size_t  i   ;
    char    tmp ;
    for (i=0; i<n; i+=4) {
        tmp     =   pt[i+3];
        pt[i+3] =   pt[i];
        pt[i]   =   tmp;

        tmp     =   pt[i+2];
        pt[i+2] =   pt[i+1];
        pt[i+1] =   tmp;
    }
}

/* ----------------------------------------------------------------------
 *  check_sac_nvhdr
 *
 *  Description: Determine the byte order of the SAC file
 *
 *  IN:
 *      const int nvhdr : nvhdr from header
 *
 *  Return:
 *      FALSE   no byte order swap is needed
 *      TRUE    byte order swap is needed
 *      -1      not in sac format ( nvhdr != SAC_HEADER_MAJOR_VERSION )
 *
  --------------------------------------------------------------------- */
static int check_sac_nvhdr(const int nvhdr)
{
    int lswap = FALSE;

    if (nvhdr != SAC_HEADER_MAJOR_VERSION) {
        byte_swap((char*) &nvhdr, SAC_DATA_SIZEOF);
        if (nvhdr == SAC_HEADER_MAJOR_VERSION)
            lswap = TRUE;
        else
            lswap = -1;
    }
    return lswap;
}

/*
 *  map_chdr_in:
 *       map strings from buffer to memory
 */
static void map_chdr_in(char *memar, char *buff)
{
    char    *ptr1;
    char    *ptr2;
    int     i;

    ptr1 = memar;
    ptr2 = buff;

    memcpy(ptr1, ptr2, 8);
    *(ptr1+8) = '\0';
    ptr1 += 9;
    ptr2 += 8;

    memcpy(ptr1, ptr2, 16);
    *(ptr1+16) = '\0';
    ptr1 += 18;
    ptr2 += 16;

    for (i=0; i<21; i++) {
        memcpy(ptr1, ptr2, 8);
        *(ptr1+8) = '\0';
        ptr1 += 9;
        ptr2 += 8;
    }
}


/* ----------------------------------------------------------------
 *  read_head_in:
 *      read sac header in and deal with possible byte swap.
 *
 *  IN:
 *      const char *name : file name, only for debug
 *      SACHEAD    *hd   : header to be filled
 *      FILE       *strm : file handler
 *
 *  Return:
 *      0   :   Succeed and no byte swap
 *      1   :   Succeed and byte swap
 *     -1   :   fail.
   -------------------------------------------------------------- */
static int read_head_in(const char *name, SACHEAD *hd, FILE *strm)
{
    char   *buffer;
    int     lswap;

    if (sizeof(float) != SAC_DATA_SIZEOF || sizeof(int) != SAC_DATA_SIZEOF) {
        fprintf(stderr, "Mismatch in size of basic data type!\n");
        return -1;
    }

    /* read numeric parts of the SAC header */
    if (fread(hd, SAC_HEADER_NUMBERS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC header %s\n", name);
        return -1;
    }

    /* Check Header Version and Endian  */
    lswap = check_sac_nvhdr(hd->nvhdr);
    if (lswap == -1) {
        fprintf(stderr, "Warning: %s not in sac format.\n", name);
        return -1;
    } else if (lswap == TRUE) {
        byte_swap((char *)hd, SAC_HEADER_NUMBERS_SIZE);
    }

    /* read string parts of the SAC header */
    if ((buffer = (char *)malloc(SAC_HEADER_STRINGS_SIZE)) == NULL) {
        fprintf(stderr, "Error in allocating memory %s\n", name);
        return -1;
    }
    if (fread(buffer, SAC_HEADER_STRINGS_SIZE, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC header %s\n", name);
        free(buffer);
        return -1;
    }
    map_chdr_in((char *)(hd)+SAC_HEADER_NUMBERS_SIZE, buffer);
    free(buffer);

    return lswap;
}

/* ============================================================== */
/* ============================================================== */
/* --------------------------------------------------------------
 *  read_sac
 *
 *  Description: Read binary SAC data from file.
 *
 *  IN:
 *      const char *name : file name
 *  OUT:
 *      SACHEAD    *hd   : SAC header to be filled
 *  Return: float pointer to the data array, NULL if failed.
 *
   ------------------------------------------------------------- */
float *read_sac(const char *name, SACHEAD *hd)
{
    FILE    *strm;
    float   *ar;
    int     lswap;
    size_t  sz;

    if ((strm = fopen(name, "rb")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", name);
        return NULL;
    }

    lswap = read_head_in(name, hd, strm);

    if (lswap == -1) {
        fclose(strm);
        return NULL;
    }

    sz = (size_t) hd->npts * SAC_DATA_SIZEOF;
    if (hd->iftype == IXY) sz *= 2;

    if ((ar = (float *)malloc(sz)) == NULL) {
        fprintf(stderr, "Error in allocating memory for reading %s\n", name);
        fclose(strm);
        return NULL;
    }

    if (fread((char*)ar, sz, 1, strm) != 1) {
        fprintf(stderr, "Error in reading SAC data %s\n", name);
        free(ar);
        fclose(strm);
        return NULL;
    }
    fclose(strm);

    if (lswap == TRUE) byte_swap((char*)ar, sz);

    return ar;
}

/* ------------------------------------------------------------
cal_dist_az_baz

DESCRIPTIONS: Calculate the great circle distance, azimuth 
              and back-azimuth in degree.

IN:        
    float evtlon: longitude of event; 
    float evtlon: longitude of event;
    float stalat: latitude of station;
    float stalon: longitude of station; 
    float res[]:  3-element array saving distance
                  and azimuth and back-aimuth. 
HISTORY:
     T. Owens, September 19, 1991
        Sept. 25 -- fixed az and baz calculations

     P. Crotwell, Setember 27, 1994
        Converted to c to fix annoying problem of fortran giving wrong
        answers if the input doesn't contain a decimal point.
     X. Feng, Jan 05, 2021
        Modified and organised the codes.
 ------------------------------------------------------------ */
void cal_dist_az_baz(float evtlat, float evtlon,\
                 float stalat, float stalon,\
                 float res[])
{
    float delta, az, baz;
    float scolat, slon, ecolat, elon;
    float a, b, c, d, e, aa, bb, cc, dd, ee, g, gg, h, hh, k, kk;
    float rhs1, rhs2, sph, rad, del, daz, dbaz, pi, piby2;

    if (fabs(stalat-evtlat)<1e-8&&fabs(stalon-evtlon)<1e-8) {
        res[0] = 0.; res[1] = 0.; res[2] = 0.;
        return;
    }

   pi    = 3.1415926535;
   piby2 = pi / 2.0;
   rad   = 2. * pi / 360.0;
/*
c scolat and ecolat are the geocentric colatitudes
c as defined by Richter (pg. 318)
c
c Earth Flattening of 1/298.257 take from Bott (pg. 3)
*/
   sph    = 1. / 298.257;
   scolat = piby2 - atan((1.-sph) * (1.-sph) * tan(stalat*rad));
   ecolat = piby2 - atan((1.-sph) * (1.-sph) * tan(evtlat*rad));
   slon   = stalon * rad;  elon = evtlon * rad;
/*
c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
c     These are defined for the pt. 1
*/
    a = sin(scolat) * cos(slon); b = sin(scolat) * sin(slon);
    c = cos(scolat); d = sin(slon); e = -cos(slon);
    g = -c * e;      h = c * d;     k = -sin(scolat);
/*
c  aa - ee are the same as a - e, except for pt. 2
*/
    aa = sin(ecolat) * cos(elon); bb = sin(ecolat) * sin(elon);
    cc = cos(ecolat); dd = sin(elon); ee = -cos(elon);
    gg = -cc * ee;    hh = cc*dd;     kk = -sin(ecolat);
/*
c  Bullen, Sec 10.2, eqn. 4
*/
    del   = acos(a*aa + b*bb + c*cc);
    delta = del/rad;
/*
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 1 is unprimed, so this is technically the baz
c
c  Calculate baz this way to avoid quadrant problems
*/
    rhs1 = (aa-d) * (aa-d) + (bb-e) * (bb-e) + cc*cc - 2.;
    rhs2 = (aa-g) * (aa-g) + (bb-h) * (bb-h) + (cc-k)*(cc-k) - 2.;
    dbaz = atan2(rhs1,rhs2);
    if (dbaz<0.0) {
        dbaz = dbaz + 2 * pi;
    }
    baz = dbaz / rad;

/*
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 2 is unprimed, so this is technically the az
*/
    rhs1 = (a-dd) * (a-dd) + (b-ee) * (b-ee) + c*c - 2.;
    rhs2 = (a-gg) * (a-gg) + (b-hh) * (b-hh) + (c-kk)*(c-kk) - 2.;
    daz  = atan2(rhs1,rhs2);
    if(daz<0.0) {
        daz = daz + 2 * pi;
    }
    az = daz / rad;

/*
c   Make sure 0.0 is always 0.0, not 360.
*/
    //if(abs(baz-360.) < .00001) baz=0.0;
    if(abs(az-360.) < .00001) az=0.0;

// Return distance in degree and azimuth in degress.
    res[0] = delta; res[1] = az; res[2] = baz;
}


/* -----------------------------------------------------------
 *
 * Delet space or line break.
 * ----------------------------------------------------------*/
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


/* ---------------------------------------------------------
 *
 * Calculate the integer power of 2 near the given integer.
 *
 * ------------------------------------------------------- */
int near_pow2(int n) {
    int m, i;
    i = (int) (log((float)n) / log(2.) + 1.);
    m = (int) pow(2., i);
    return m;
}



/* -------------------------------------------------------
 *
 * Cross-psectral beam-forming with given frequency.
 *
 * ---------------------------------------------------- */
void array_spectra( ARRAYPARA p, fftw_complex **sp, char *system) {
    int i, j, k, m, nccf;
    float tmp1, tmp2, slow, baz, ds, db, shift;
    float **rpos;
    fftw_complex *ccfd;

    nccf = p.nsac * (p.nsac-1) / 2;
    ccfd = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * nccf );
    rpos = (float **) malloc( sizeof(float *) * nccf );
    for(i = 0; i < nccf; i ++)
        rpos[i] = (float *) malloc( sizeof(float) * 2 );

// Obtain relative coordinates of seismic array.
    if( !(strcmp("lonlat", system))) {
        const float g = 111.19492664455;
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

// Initialization of array processing spectra.
    for(i = 0; i < p.ns; i ++)
        for(j = 0; j < p.nb; j ++) {
            sp[i][j][0] = 0.;
            sp[i][j][1] = 0.;
        }

// Scanning slowness and back-azimuth over given frequency band.
    float df, f;
    df = (p.f2 - p.f1) / (p.fn-1);
    for(m = 0; m < p.fn; m ++) {
        f = p.f1 + df * m;
        printf("Scanning at frequency %.5f Hz ...\n", f);
        tmp1 = -2.* PI * f;
        ds   = (p.s2 - p.s1) / (p.ns - 1);
        db   = (p.b2 - p.b1) / (p.nb - 1);
    // Compute cross-spectra.
        k = 0;
        for(i = 0; i < p.nsac-1; i ++)
            for(j = i+1; j < p.nsac; j ++){
                ccfd[k][0] = p.fft[m][i][0] * p.fft[m][j][0] \
                           + p.fft[m][i][1] * p.fft[m][j][1]; 
                ccfd[k][1] = p.fft[m][j][0] * p.fft[m][i][1] \
                           - p.fft[m][i][0] * p.fft[m][j][1];
                k ++;
        }
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
    }

// Release dynamic memories.
    for(i = 0; i < nccf; i ++) {
        free(rpos[i]);
    }
    free(rpos); fftw_free(ccfd);
}
