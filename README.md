## Seismic data
    Five data sets are prepared in folder `Data`.
### 1 DATA1
    Synthetic vertical-component seismograms with YASEIS (Ma, 2013)
        based on the PREM.
    The source mechanism is explosion at depth of 20 km and the
        average distance is about 20 degress.

### 2 DATA2
    Synthetic vertical-component seismograms with YASEIS (Ma, 2013)
        based on the PREM.
    The source mechanism is explosion at depth of 20 km and the
        average distance is about 50 degress.

### 3 DATA3
    Synthetic vertical-component seismograms with CPS (Herrmann, 2013)
        based on a multi-layered model.
    The source mechanism is explosion at depth of 36 km and the
        average distance is about 500 km.

### 4 DATA4
    Real vertical-component seismograms of NO network from NORSAR data
        center in Norway.
    The major array (NORES) is composed of 7 sub-arrays.

### 5 DATA5
    Real vertical-component seismograms of NO network from NORSAR data
        center in Norway.
    The circle-shaped array is named ARCES and composed of 25 stations.


## Example
### Input file
    The contents of input file are described by the following lines.
`
listfile1    #File of SAC list names;
200.         #Lower limit of time in second;
400.         #Upper limit of time in second;
0.25         #Lower limit of frequency of plane wave;
0.26         #Upper limit of frequency of plane wave;
0.01         #Lower limit of slowness in s/km;
0.15         #Upper limit of slowness in s/km;
0.           #Lower limit of back-azimuth in degree;
360.         #Upper limit of back-azimuth in degree;
81           #Segments of slowness;
181          #Segments of bazk-azimuth;
lonlat       #Coordinate system: "lonlat", "km" or "meter";
out.txt      #File of results.
`
### Steps
    list SAC format file in a file, here `listfile1`;
    set parameters in input file, here `input.para1`;
    run an example with `../bin/cross_spec_beam input.para1`;
    plot the spectra with `python plot_spectra.py input.para1 deg`.

## Installztion
    Go into folder `src` and run `make` to compile and executable
    command `cross_spec_beam` will be generated in folder `bin`.

## Dependencies
    `gcc` or `icc` compiler;
    `fftw3` for doing FFT;
    `python3` for visualization.
