module Transforms
using ..WT
export cwtft
function cwtft{T <: Real}(Y :: AbstractArray{T}, dt :: Number;
			pad :: Bool = false, dj :: Number = 0.25, s0 :: Number = 2.*dt, J1 :: Number = -1,
			mother = WT.morlet, param :: Number = WT.sparam(mother),
			freq :: Array{Float64,1} = [])

		#Y=Y[:];
		n1 = length(Y);

		if J1 == -1
		        J1=floor(Int,(log(n1*dt/s0)/log(2.))/dj);
		end
		#....construct time series to analyze, pad if necessary
		x = Y - mean(Y);
		if pad
		        base2 = floor(Int,log(n1)/log(2) + 0.4999);   # power of 2 nearest to N
		        x = [x, zeros(2^(floor(Int,base2)+1)-n1)];
		end
		n = length(x);

		#....construct wavenumber array used in transform [Eqn(5)]
		k = 1:floor(Int,n/2);
		k = [0.;  k;  -k[floor(Int,(n-1)/2):-1:1]]*((2*pi)/(n*dt));
		#....compute FFT of the (padded) time series
		f = fft(x);    # [Eqn(3)]
		#....construct SCALE array & empty PERIOD & WAVE arrays
		if isempty(freq)
			scale = s0*2.^((0:J1)*dj);
		else
			scale = 1./(FourierFactor(mother,param)*freq);
			J1 = length(scale)-1;
		end

		wave = zeros(Complex{Float64},J1+1,n);  # define the wavelet array


		# loop through all scales and compute transform
		for a1 in 1:J1+1
		        daughter = WT.Daughter(mother,scale[a1],k,param,n)
		        wave[a1,:] = ifft(f.*daughter)  # wavelet transform[Eqn(4)]
		end
		fourier_factor=WT.FourierFactor(mother,param);
		period = fourier_factor*scale;
		coi = WT.COI(mother,fourier_factor).*dt*[1E-5; 1:((n1+1)/2-1); flipdim((1:(n1/2-1)),1); 1E-5];  # COI [Sec.3g]
		wave = wave[:,1:n1];  # get rid of padding before returning


		return wave,period,scale,coi
end
end
