## F-K cross-spectral beamforming
    Searching the back-azimuth and corresponding apparent slowness (or ray parameter)    
    with F-K cross-spectral beamforming array method.
## Seismic data tested
    Several data sets are prepared in folder "Data".
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
    
### 6 DATA_EVENT
    Earthquake event data: 2017-06-22T12:31:04|13.7527|-90.9488|46.82|at,pt,us|NEIC PDE|us|pt17173000,at00ory83s,us20009p1a|mww|6.8|us|NEAR COAST OF GUATEMALA
    


### Usage tips
#### input file preparation
The contents of input file are described by the following lines.

`listfile1    #File of SAC list names;`    
`200.         #Lower limit of time in second;`    
`400.         #Upper limit of time in second;`    
`0.25         #Lower limit of frequency of plane wave;`    
`0.26         #Upper limit of frequency of plane wave;`    
`0.01         #Lower limit of slowness in s/km;`    
`0.15         #Upper limit of slowness in s/km;`    
`0.           #Lower limit of back-azimuth in degree;`    
`360.         #Upper limit of back-azimuth in degree;`    
`81           #Segments of slowness;`    
`181          #Segments of bazk-azimuth;`    
`lonlat       #Coordinate system: "lonlat", "km" or "meter";`    
`out.txt      #File of results.`

#### steps
list SAC format files in a file, here 
`listfile1`;
set parameters in input file, here 
`input.para1`;
run an example with    
`../bin/cross_spec_beam input.para1`    
plot the spectra with      
`python plot_spectra.py input.para1 deg`.


***
## Back-azimuth and ray parameter of P or pP phase
`2017-06-22T12:31:04|13.7527|-90.9488|46.82|at,pt,us|NEIC PDE|us|pt17173000,at00ory83s,us20009p1a|mww|6.8|us|NEAR COAST OF GUATEMALA`     

#### waveforms
![wave](https://github.com/geophydog/Beamforming_in_frequency_domain/blob/main/Data/images/wave.png)

#### beamfoming power
`../bin/cross_spec_beam input.para9`       
`python plot_spectra.py input.para9 deg`        
![beam](https://github.com/geophydog/Beamforming_in_frequency_domain/blob/main/Data/images/630-700.png)


## Installation
Go into the directory `src` and type `make` to compile and the executable
command `cross_spec_beam` will be generated in the upper-level directory `bin`.

## Dependencies
`gcc` or `icc` compiler;
`fftw3` for doing FFT;
`python3` for visualization.
