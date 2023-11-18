"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""


import OrbitalElements
using Printf
using HDF5

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64    = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)


# params to specify
rmin,rmax,rc,da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε = 0.,Inf,1.0,0.001,0.001,0.001,32,0.01,0.001,100,1.e-6,1.e6,1.e-10

# these can be repurposed if one wants to run convergence tests with NINT or EDGE.
NINTmin = 32
NINTstep = 10
NINTsteps = 1#00

EDGEmin = -4.0
EDGEmax = log10(0.33)
EDGEsteps = 1#00
EDGEstep = (EDGEmax-EDGEmin)/(EDGEsteps)

Ω₁accuracy = zeros(Float64,NINTsteps,EDGEsteps)
Ω₂accuracy = zeros(Float64,NINTsteps,EDGEsteps)
Ω₁e = zeros(Float64,NINTsteps,EDGEsteps)
Ω₂e = zeros(Float64,NINTsteps,EDGEsteps)

Usteps = 200
dU = 2/(Usteps-1)
Θgrid = zeros(Float64,NINTsteps,EDGEsteps,Usteps)

Uvals = zeros(Float64,Usteps)
for Ui=1:Usteps
    Uvals[Ui] = -1 + (Ui-1)*dU
end



# select an (a,e) value for the orbit
a,e = 0.001, 0.99

for a in [0.001,1.0,1000.0]
    for e in [0.01,0.5,0.99]
        Ngrid = zeros(Float64,NINTsteps,EDGEsteps)
        Egrid = zeros(Float64,NINTsteps,EDGEsteps)

        # OrbitalElements.ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)
        # number of steps inside EDGE is calculable

        # now let's loop through and see how accurate we can get
        for NINTi=1:NINTsteps
            nval = (NINTi-1)*NINTstep+NINTmin
            for EDGEi=1:EDGEsteps
                eval = 10^((EDGEi-1)*EDGEstep+EDGEmin)
                Ngrid[NINTi,EDGEi] = nval
                Egrid[NINTi,EDGEi] = eval
                params = OrbitalElements.OrbitalParameters(Ω₀,rmin,rmax,rc,eval,TOLECC,TOLA,nval,da,de,ITERMAX,invε)

                Ω₁accuracy[NINTi,EDGEi],Ω₂accuracy[NINTi,EDGEi] = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,params)
                Ω₁e[NINTi,EDGEi],Ω₂e[NINTi,EDGEi] = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)

                for Ui=1:Usteps
                    uval = -1 + (Ui-1)*dU
                    Θgrid[NINTi,EDGEi,Ui] = OrbitalElements.ΘAE(ψ,dψ,d2ψ,uval,a,e,params)
                end
                #println("oldae O1={$Ω₁accuracy[NINTi,EDGEi]} O2=$Ω₂r")
            end
        end


        h5open("../figureA1/ValidationDataA$(a)E$(e).h5", "w") do file
            # Mappings parameters
            write(file, "o1accuracy",Ω₁accuracy)
            write(file, "o2accuracy",Ω₂accuracy)
            write(file, "ngrid",Ngrid)
            write(file, "egrid",Egrid)
            write(file, "o1true",Ω₁e)
            write(file, "o2true",Ω₂e)
            write(file, "tgrid",Θgrid)
            write(file, "uvals",Uvals)
        end
    end # end e loop
end # end a loop
