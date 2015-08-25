function wavelet(Y,dt; pad=0,dj=0.25,s0=2*dt,J1=-1,mother="MORLET",param=-1,parallel=false)

Y=Y[:];
n1 = length(Y);

if J1 == -1
	 J1=int((log(n1*dt/s0)/log(2))/dj);
end
#....construct time series to analyze, pad if necessary
x = Y - mean(Y);
if pad == 1
	base2 = int(log(n1)/log(2) + 0.4999);   # power of 2 nearest to N
	x = [x, zeros(2^(int(base2)+1)-n1)];
end
n = length(x);

#....construct wavenumber array used in transform [Eqn(5)]
k = 1:int(n/2);
k = k*((2*pi)/(n*dt));
k = [0.;  k;  -k[(int((n-1)/2)-1):-1:1]];
#....compute FFT of the (padded) time series
# ft=plan_fft(x);
f = fft(x);    # [Eqn(3)]
# ift=plan_ifft(f);
#....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);
period = scale;
fourier_factor=0.;
wave = zeros(Complex,J1+1,n);  # define the wavelet array
    
# loop through all scales and compute transform
if parallel
        
        daughter, fourier_factor, coi, dofmin=wave_bases(mother,k,scale[1],param);
        wave[1,:]=ifft(f.*daughter);
        waver=SharedArray(Float64,(J1,n), init= waver -> 0);
        waveim=SharedArray(Float64,(J1,n), init= waver -> 0);
        @sync @parallel for a1 in 1:J1
            daughter, ff, coi, dofmin=wave_bases(mother,k,scale[a1+1],param);
            temp=ifft(f.*daughter); # wavelet transform[Eqn(4)]
            waver[a1,:]=real(temp); 
            waveim[a1,:]=imag(temp); 
        end
        wave[2:end,:]=complex(sdata(waver),sdata(waveim));
else
    coi=zeros(size(x));
    for a1 in 1:J1+1
        daughter, fourier_factor, coi, dofmin=wave_bases(mother,k,scale[a1],param)
        wave[a1,:] = ifft(f.*daughter);  # wavelet transform[Eqn(4)]
    end
end


#     println( typeof(fourier_factor))
#     println( fourier_factor)
period = fourier_factor*scale;
coi = coi.*dt*[1E-5; 1:((n1+1)/2-1); flipdim((1:(n1/2-1)),1); 1E-5];  # COI [Sec.3g]
wave = wave[:,1:n1];  # get rid of padding before returning


	return wave,period,scale,coi
end
