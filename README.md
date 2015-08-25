# wave_julia
Continuous Wavelet Transform in Julia based on Torrence and Compo (1998) Script, which is available on: http://paos.colorado.edu/research/wavelets/
## Use

To compute the Wavelet transform of a signal `x` with temporal spacing `dt` simply type:

```julia
wave,period,scales,coi=wavelet(x,dt)
```

which gives you:

The wavelet coefficient matrix `wave` (dimensions: `length(scales) x length(x)`).

The fourier periods vector `period` (dimensions: `length(scales)`).

The scales vector `scales` containing the optimal scales for the wavelet transform after Torrence and Compo (1998) [eq. 9].

The cone of influence vector `coi` (dimensions: `length(x)`) to identify coefficients under influence of border effects.

##Optional input (with default values):

`pad=false` : enable zero padding for fft speed increase

`dj=0.25` : scale spacing coeff.

`s0=2*dt` : smallest scale

`J1=-1` : largest scale, if `-1` the largest scale is computed after Torrence and Compo (1998) [eq. 10]

`mother="MORLET"` : Wavelet function, possible Wavelets are: `MORLET`, `PAUL`, `DOG`

`param=-1` : Parameter for Wavelet function 

`parallel=false` : see next section

For more details read Torrence and Compo (1998).


## Parallel
The code is parallel for a single node using SharedArrays, which is significantly slower for small size time series (10^3 data points), but faster for large time series (larger than 10^4). Simply use argument keyword: 

```julia
parallel=true
```

#References

C. Torrence und G. Compo. A practical guide to wavelet analysis. *Bulletin of the Ame-
rican Meteorological Society*, 79:61–78, 1998. doi: 10.1175/1520-0477(1998)079 0061:
APGTWA 2.0.CO;2.


