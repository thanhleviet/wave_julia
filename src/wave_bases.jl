function wave_bases(mother,k,scale,param)

mother = uppercase(mother);
n = length(k);

if mother == "MORLET"  #-----------------------------------  Morlet
	if param == -1
		 param = 6;
	end
	k0 = param;
    expnt = -(scale*k - k0).^2/2;
		# expnt[k.<0]=0.;
	# for i in 1:n
	# 	if k[i] <= 0
	# 		expnt[i] = 0
	# 	end
	# end
	norm = sqrt(scale*k[2])*(pi^(-1//4))*sqrt(n);    # total energy=N   [Eqn(7)]
	daughter = norm*exp(expnt);
	daughter[k.<0]=0.;
	# for i in 1:n
	# 	if k[i] <= 0
	# 		daughter[i] = 0
	# 	end
	# end
	fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); # Scale-->Fourier [Sec.3h]
	coi = fourier_factor/sqrt(2);                  # Cone-of-influence [Sec.3g]
	dofmin = 2;                                    # Degrees of freedom
elseif mother == "PAUL"  #--------------------------------  Paul
	if param == -1
		param = 4;
	end
	m = param;
	expnt = -(scale*k);
	# for i in 1:n
	# 	if k[i] <= 0
	# 		expnt[i] = 0
	# 	end
	# end
	norm = sqrt(scale*k[2])*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
	daughter = norm*((scale*k).^m).*exp(expnt);
	# for i in 1:n
	# 	if k[i] <= 0
	# 		daughter[i] = 0
	# 	end
	# end
	daughter[k.<0]=0.;
	fourier_factor = 4*pi/(2*m+1);
	coi = fourier_factor*sqrt(2);
	dofmin = 2;
elseif mother == "DOG"  #--------------------------------  DOG
	if param == -1
		param = 2;
	end
	m = param;
	expnt = -(scale*k).^2/2.0;
	# for i in 1:n
	# 	if k[i] <= 0
	# 		expnt[i] = 0
	# 	end
	# end
	norm = sqrt(scale*k[2]/gamma(m+0.5))*sqrt(n);
	daughter = -norm*(i^m)*((scale*k2)^m).*exp(expnt);
	# for i in 1:n
	# 	if k[i] <= 0
	# 		daughter[i] = 0
	# 	end
	# end
	daughter[k.<0]=0.;
	fourier_factor = 2*pi*sqrt(2./(2*m+1));
	coi = fourier_factor/sqrt(2);
	dofmin = 1;
else
	error("Mother must be one of MORLET,PAUL,DOG")
end

return daughter,fourier_factor,coi,dofmin

end
