# wave_julia
Continuous Wavelet Transform in Julia based on Torrence and Compo (1998) Script, which is available on: http://paos.colorado.edu/research/wavelets/
## Use

To compute the Wavelet transform (using "optimum" set of scales) of a signal `x` with temporal spacing `dt` simply type:
```julia
using wave_julia
```

```julia
wave,period,scales,coi=wave_julia.Transforms.cwtft(x,dt)
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


For more details read Torrence and Compo (1998).




##References

C. Torrence and G. Compo. A practical guide to wavelet analysis. *Bulletin of the American Meteorological Society*, 79:61–78, 1998. doi: 10.1175/1520-0477(1998)079 0061:APGTWA 2.0.CO;2.


