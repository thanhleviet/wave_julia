# wave_julia
Continuous Wavelet Transform in Julia based on Torrence and Compo (1998) http://paos.colorado.edu/research/wavelets/

Parallel for single node using SharedArrays, which is significantly \slower\ for small size time series (10^3 data points), but faster for large time series (larger than 10^4). Simply use argument keyword: 

parallel=true

