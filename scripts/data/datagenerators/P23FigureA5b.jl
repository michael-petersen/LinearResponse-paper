
using HDF5

using FiniteHilbertTransform
using OrbitalElements

const bc, M, G = 1.,1. ,1.

ψ(r::Float64)::Float64       = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64    = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψPlummer(r,bc,M,G)
Ω₀      =    OrbitalElements.Ω₀Plummer(bc,M,G)


Ku = 200
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

n1 = -1
n2 = 2

# make a blank array of values to test with omg (map the whole space?)
omgsamples = -10 .^collect(-4:0.001:-1)
nsamples = size(omgsamples)[1]
Darray = zeros(Complex{Float64}, nsamples,Ku)
Parray = zeros(Complex{Float64}, nsamples,Ku)
Qarray = zeros(Complex{Float64}, nsamples,Ku)

romg = 0.0
restag = "-12"
model = "plummer"

for (i,omg) in enumerate(omgsamples)
    omgin = romg+omg*1im
    varpi = OrbitalElements.Getϖ(omgin/Ω₀,n1,n2,dψ,d2ψ) # Getting the rescaled frequency from # Dimensionless frequency rescaled by Omega0
    FiniteHilbertTransform.GettabD!(varpi,FHT)
    Darray[i,:] = FHT.tabDLeg
    Parray[i,:] = FHT.tabPLeg
    Qarray[i,:] = FHT.tabQLeg
end


h5open("../figureA5/DConvergence-$model-K$Ku-$restag-$romg.h5", "w") do file
    write(file,"tabDLeg",Darray)
    write(file,"tabPLeg",Parray)
    write(file,"tabQLeg",Qarray)
end
