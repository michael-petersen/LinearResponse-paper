
import OrbitalElements
using Printf
using HDF5

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.

ψ(r::Float64)::Float64    = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

# make a HIGH RES version of the frequencies

Amin = -3
Amax = 3
Asteps = 1000
Astep = (Amax-Amin)/(Asteps-1)

Emin = 0.001
Emax = 0.999
Esteps = 1000
Estep = (Emax-Emin)/(Esteps-1)


function run_inversion_grid(Amin,Amax,Astep,Emin,Emax,Estep)
    # params to specify
    rmin,rmax,rc,da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε = 0.,Inf,1.0,0.001,0.001,0.001,32,0.01,0.001,100,1.e-6,1.e6,1.e-10

    # now we need to make a grid of re-inversions to check accuracy there
    Ω₁inversion = zeros(Float64,Asteps,Esteps)
    Ω₂inversion = zeros(Float64,Asteps,Esteps)
    αe = zeros(Float64,Asteps,Esteps)
    βe = zeros(Float64,Asteps,Esteps)
    Agrid = zeros(Float64,Asteps,Esteps)
    Egrid = zeros(Float64,Asteps,Esteps)


    for Ai=1:Asteps
        aval = 10^((Ai-1)*Astep+Amin)
        for Ei=1:Esteps
            eval = (Ei-1)*Estep+Emin
            Agrid[Ai,Ei] = aval
            Egrid[Ai,Ei] = eval
            params = OrbitalElements.OrbitalParameters(Ω₀,rmin,rmax,rc,EDGE,TOLECC,TOLA,NINT,da,de,ITERMAX,invε)
            αe[Ai,Ei],βe[Ai,Ei] = OrbitalElements.αβFromAE(ψ,dψ,d2ψ,aval,eval,params)
            Ω₁inversion[Ai,Ei],Ω₂inversion[Ai,Ei] = OrbitalElements.ComputeAEFromαβ(ψ,dψ,d2ψ,αe[Ai,Ei],βe[Ai,Ei],params)
        end
    end


    h5open("../figureA3/InversionData_PlummerRelative.h5", "w") do file
        write(file, "afreq",αe)
        write(file, "bfreq",βe)
        write(file, "agrid",Agrid)
        write(file, "egrid",Egrid)
        write(file, "aest",Ω₁inversion)
        write(file, "eest",Ω₂inversion)
    end
end

@time run_inversion_grid(Amin,Amax,Astep,Emin,Emax,Estep)
