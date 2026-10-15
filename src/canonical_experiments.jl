# Small, explicitly identified models for the canonical experiment collection.

# These are numerical experiments, not replications of all cited devices. The
# membrane is a prestressed axisymmetric membrane reduction, the Helmholtz tube
# is a passive transmission-line network, and the prism is a causal effective
# medium. Each result records that distinction. No Pluto or rendering dependency.
export CanonicalRecord, membrane_model, membrane_response, membrane_record,
       helmholtz_record, helmholtz_modulus, canonical_spatial_record

using SparseArrays

struct CanonicalRecord
    x::Vector{Float64}
    y::Vector{Float64}
    times::Vector{Float64}
    field::Matrix{Float32} # x fastest, then y, then measurement time
    valid::BitVector
    internal::Matrix{Float32}
    diagnostics::Dict{String,Float64}
end

_burst(t,f,cycles=6.0)=0<t<cycles/f ? sin(2π*f*t)*sinpi(t*f/cycles)^2 : 0.0

"""Axisymmetric, tensioned annulus with a rigid central disk. SI units.
The tension is our specimen choice, NOT an inferred parameter of Yang 2008.
"""
function membrane_model(;radius=.01,disk_radius=.003,added_mass=.0003,
                        surface_density=.2744,tension=5.0,elements=32)
    0<disk_radius<radius && elements>=4 && added_mass>=0 &&
        surface_density>0 && tension>0 || throw(ArgumentError("invalid membrane"))
    r=collect(range(disk_radius,radius;length=elements+1)); n=elements
    M=zeros(n,n); K=zeros(n,n); b=zeros(n)
    # Two-point Gauss exactly integrates the radial consistent mass and load.
    for e in 1:elements
        h=r[e+1]-r[e]
        for ξ in (-inv(sqrt(3)),inv(sqrt(3)))
            N=[(1-ξ)/2,(1+ξ)/2]; dN=[-1/h,1/h]
            weight=π*h*(N[1]*r[e]+N[2]*r[e+1])
            for i in 1:2
                a=e+i-1; a>n && continue
                b[a]+=weight*N[i]
                for j in 1:2
                    c=e+j-1; c>n && continue
                    M[a,c]+=weight*surface_density*N[i]*N[j]
                    K[a,c]+=weight*tension*dN[i]*dN[j]
                end
            end
        end
    end
    disk_area=π*disk_radius^2; area=π*radius^2
    M[1,1]+=added_mass+surface_density*disk_area; b[1]+=disk_area
    (;M,K,b,r,area,radius,disk_radius,added_mass,tension,surface_density)
end

function membrane_response(m,f;back_depth=0.0,loss_rate=0.0)
    f>0 && loss_rate>=0 && back_depth>=0 || throw(ArgumentError("invalid frequency/loss"))
    ρ=1.2; c=343.; Z=ρ*c; ω=2π*f
    K=m.K+(back_depth>0 ? (ρ*c^2/(m.area*back_depth))*(m.b*m.b') : zero(m.K))
    C=loss_rate*m.M
    D=ComplexF64.(K-ω^2*m.M)-im*ω*C
    compliance=dot(m.b,D\m.b)/m.area
    mass=-1/(ω^2*compliance)
    impedance=inv(-im*ω*compliance)
    if back_depth>0
        r=(impedance-Z)/(impedance+Z); t=0im
    else
        r=impedance/(impedance+2Z); t=2Z/(impedance+2Z)
    end
    (;mass,impedance,r,t,absorption=1-abs2(r)-abs2(t))
end

function _membrane_history(m,f;back_depth=0.0,loss_rate=.25,duration=16/f,steps_per_cycle=100)
    ρ=1.2; c=343.; Z=ρ*c; ports=back_depth>0 ? 1 : 2
    K=m.K+(back_depth>0 ? (ρ*c^2/(m.area*back_depth))*(m.b*m.b') : zero(m.K))
    C=loss_rate*m.M+(ports*Z/m.area)*(m.b*m.b')
    dt=1/(steps_per_cycle*f); times=collect(0:dt:duration); n=length(m.b)
    u=zeros(n); v=zeros(n); a=zeros(n); U=zeros(Float32,n,length(times)); V=zeros(length(times))
    factor=cholesky(Symmetric(m.M+dt/2*C+dt^2/4*K))
    work=0.; diss=0.; lastforce=zeros(n)
    for j in 2:length(times)
        up=u+dt*v+dt^2/4*a; vp=v+dt/2*a
        force=2*_burst(times[j],f)*m.b
        an=factor\(force-C*vp-K*up)
        un=up+dt^2/4*an; vn=vp+dt/2*an
        midpoint=(vn+v)/2
        work+=dot((force+lastforce)/2,midpoint)*dt
        diss+=dot(midpoint,C*midpoint)*dt
        u=un; v=vn; a=an; lastforce=force
        U[:,j]=u; V[j]=dot(m.b,v)/m.area
    end
    energy=(dot(v,m.M*v)+dot(u,K*u))/2
    balance=abs(energy+diss-work)/max(abs(work),eps())
    (;times,U,V,energy,work,diss,balance,dt)
end

"""True Newmark transient with exact plane-wave port radiation.
Field in the ducts is evaluated by retarded port signals, not phase rotation.
"""
function membrane_record(m,f;back_depth=0.0,loss_rate=.25,wall=false)
    history=if wall
        dt=1/(100*f); times=collect(0:dt:16/f)
        (times=times,U=zeros(Float32,length(m.b),length(times)),V=zeros(length(times)),
            energy=0.,work=0.,diss=0.,balance=0.,dt=dt)
    else
        _membrane_history(m,f;back_depth,loss_rate)
    end
    t=history.times; dt=history.dt; c=343.; Z=1.2*c
    x=collect(range(-1.25,1.25;length=101)); y=[0.0]
    field=zeros(Float32,length(x),length(t)); valid=trues(length(x))
    delayed(a,τ)=τ<0 || τ>last(t) ? 0.0 : begin
        q=τ/dt+1; lo=clamp(floor(Int,q),1,length(t)); hi=min(lo+1,length(t))
        a[lo]+(q-lo)*(a[hi]-a[lo])
    end
    reflected=[_burst(τ,f)-(wall ? 0 : Z*history.V[j]) for (j,τ) in enumerate(t)]
    transmitted=wall || back_depth>0 ? zeros(length(t)) : Z*history.V
    for j in eachindex(t),i in eachindex(x)
        field[i,j]=x[i]<0 ? _burst(t[j]-x[i]/c,f)+delayed(reflected,t[j]+x[i]/c) :
            delayed(transmitted,t[j]-x[i]/c)
    end
    back_depth>0 && (valid[x.>0].=false)
    response=membrane_response(m,f;back_depth,loss_rate)
    diag=Dict("frequency_Hz"=>f,"time_step_s"=>dt,"surface_mass_real_kg_m2"=>real(response.mass),
        "reflectance"=>wall ? 1. : abs2(response.r),"transmittance"=>wall ? 0. : abs2(response.t),
        "absorption"=>wall ? 0. : response.absorption,"energy_balance_relative"=>history.balance,
        "source_work_J"=>history.work,"membrane_tension_N_m"=>m.tension)
    wall && delete!(diag,"surface_mass_real_kg_m2") # rigid control has zero acceleration
    CanonicalRecord(x,y,t,field,valid,wall ? zero(history.U) : history.U,diag)
end

"""Dynamic bulk modulus of a lossless-fluid, damped-neck network, e^-iωt."""
function helmholtz_modulus(f;open=true,ρ=1000.,c=1480.,area=16e-6,pitch=.0092,
                          cavity=62.8e-9,neck_length=.001,neck_area=π*.0005^2,loss=.015)
    open || return complex(ρ*c^2)
    L=ρ*neck_length/neck_area; C=cavity/(ρ*c^2); ω=2π*f
    β=1/(ρ*c^2)+C/(area*pitch)/(1-ω^2*L*C-im*loss*ω*sqrt(L*C))
    inv(β)
end

"""Passive pressure/volume-flow transmission line with real Helmholtz DOFs.
Implicit midpoint has an exact discrete work/dissipation identity.
"""
function helmholtz_record(f;open=true,cells=80,loss=.015)
    cells>=40 && f>0 && loss>=0 || throw(ArgumentError("invalid duct"))
    ρ=1000.; c=1480.; area=16e-6; dx=.0092; B=ρ*c^2
    C=area*dx/B; L=ρ*dx/area; Z=ρ*c/area
    branches=open ? collect(30:41) : Int[]; nb=length(branches); n=2cells-1+2nb
    mass=vcat(fill(C,cells),fill(L,cells-1),fill(ρ*.001/(π*.0005^2),nb),fill(62.8e-9/B,nb))
    A=spzeros(n,n)
    for i in 1:cells-1
        q=cells+i; A[i,q]=-1; A[i+1,q]=1; A[q,i]=1; A[q,i+1]=-1
    end
    A[1,1]=-1/Z; A[cells,cells]=-1/Z
    Lb=ρ*.001/(π*.0005^2); Cb=62.8e-9/B; R=loss*sqrt(Lb/Cb)
    for (j,i) in enumerate(branches)
        q=2cells-1+j; p=2cells-1+nb+j
        A[i,q]=-1; A[q,i]=1; A[q,p]=-1; A[p,q]=1; A[q,q]=-R
    end
    dt=1/(100*f); times=collect(0:dt:1.2e-3); D=spdiagm(0=>mass)
    left=lu(D-dt/2*A); right=D+dt/2*A
    state=zeros(n); field=zeros(Float32,cells,length(times)); internal=zeros(Float32,nb,length(times))
    work=0.; dissipated=0.
    for j in 2:length(times)
        forcing=zeros(n); forcing[1]=(_burst(times[j],f)+_burst(times[j-1],f))/Z
        next=left\(right*state+dt*forcing); mid=(state+next)/2
        work+=dt*dot(forcing,mid)
        dissipated+=dt*(mid[1]^2/Z+mid[cells]^2/Z+R*sum(abs2,view(mid,2cells:2cells-1+nb)))
        state=next; field[:,j]=view(state,1:cells)
        nb>0 && (internal[:,j]=view(state,2cells+nb:n))
    end
    energy=dot(state,mass.*state)/2
    diag=Dict("frequency_Hz"=>f,"time_step_s"=>dt,
        "bulk_modulus_real_Pa"=>real(helmholtz_modulus(f;open,loss)),
        "energy_balance_relative"=>abs(energy+dissipated-work)/max(abs(work),eps()),
        "source_work_J"=>work,"resonance_Hz"=>1/(2π*sqrt(Lb*Cb)))
    CanonicalRecord(collect(((1:cells).-0.5).*dx),[0.0],times,field,trues(cells),internal,diag)
end

include("canonical_spatial.jl")
