module WT
export name,sparam,FourierFactor,COI,Morlet,Paul,DOG

abstract ContinuousWavelet

for (TYPE, NAMEBASE, PARAM, FFFUNCTION, COIFUNCTION, DFUNCTION, SUPERCLASS) in (
    (:Morlet, "morlet",
            6., #PARAM
            :((4*π)/(m + sqrt(2 + m^2))), #FourierFactorFunction
            :(1/sqrt(2)), #COIFunction
            :(daughter=sqrt(scale*k[2])*(π^(-1//4))*sqrt(n)*exp(-(scale*k - m).^2/2)), #DaughterFunction
            :ContinuousWavelet),
    (:Paul, "paul",
            4.,  #PARAM
            :(4*pi/(2*m+1)), #FourierFactorFunction
            :(1*sqrt(2)), #COIFunction
            :(daughter=sqrt(scale*k[2])*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n)*exp(-(scale*k))), #DaughterFunction
            :ContinuousWavelet),
    (:DOG, "dog",
            2.,  #PARAM
            :(2*pi*sqrt(2./(2*m+1))), #FourierFactorFunction
            :(1/sqrt(2)), #COIFunction
            :(expnt = -(scale*k).^2/2.0;
              norm = sqrt(scale*k[2]/gamma(m+0.5))*sqrt(n);
              daughter = -norm*(i^m)*((scale*k2)^m).*exp(expnt)), #DaughterFunction
            :ContinuousWavelet) # TODO moments
        )
    @eval begin
        immutable $TYPE <: $SUPERCLASS end
        name(::$TYPE) = string($NAMEBASE)::ASCIIString
        sparam(::$TYPE) = $PARAM
        FourierFactor(::$TYPE,m::Real)= $FFFUNCTION
        COI(::$TYPE,FFactor::Real)=FFactor* $COIFUNCTION
        function Daughter(::$TYPE,scale::Real,k::Array{Float64,1},m::Real,n::Integer)
            $DFUNCTION
           daughter[k.<0]=0.;
           return daughter
        end
    end
    CONSTNAME = symbol(NAMEBASE)
    @eval begin
        const $CONSTNAME = $TYPE()                  # type shortcut
    end
end

end
