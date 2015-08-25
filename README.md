# wave_julia
Continuous Wavelet Transform in Julia based on Torrence and Compo (1998) Script, which is available on: http://paos.colorado.edu/research/wavelets/

## Parallel
The code is parallel for a single node using SharedArrays, which is significantly \slower\ for small size time series (10^3 data points), but faster for large time series (larger than 10^4). Simply use argument keyword: 

```julia
parallel=true
```


