# Shared P2 displacement / P1 volumetric-stress forms for concentric solids.
# q is independent on the two material sides, displacement is perfectly bonded.
function _radial_elastic_forms(s,n,nodes,elements,materials)
    nv=length(nodes); components=n==0 ? 1 : 2; nu=components*nv
    np=length(elements)+1+(s.scenario==:coated ? 1 : 0)
    I=Int[]; J=Int[]; mass=Float64[]; dev=Float64[]; coat=Float64[]
    Ib=Int[]; Jb=Int[]; bv=Float64[]; Iq=Int[]; Jq=Int[]; qv=Float64[]
    mat=_coated_materials(s)
    for (e,ids) in enumerate(elements)
        ρ,K,μ=getproperty(mat,materials[e]); K=real(K); μ=real(μ)
        lo,hi=nodes[ids[1]],nodes[ids[3]]; h=hi-lo
        shift=e>s.steel_elements ? 1 : 0; qids=[e+shift,e+1+shift]
        uids=components==1 ? collect(ids) : vcat(collect(ids),nv .+ collect(ids))
        ne=length(uids); me=zeros(ne,ne); ke=zeros(ne,ne); be=zeros(ne,2); qe=zeros(2,2)
        for (ξ,w) in zip(_RADIAL_GAUSS_X,_RADIAL_GAUSS_W)
            r=(lo+hi)/2+h*ξ/2; weight=w*h*r/2
            N=[ξ*(ξ-1)/2,1-ξ^2,ξ*(ξ+1)/2]; dN=[ξ-.5,-2ξ,ξ+.5]*2/h
            Q=[(1-ξ)/2,(1+ξ)/2]; B=zeros(4,ne)
            B[1,1:3]=dN; B[2,1:3]=N/r; B[3,1:3]=-n*N/r
            if components==2
                B[2,4:6]=n*N/r; B[3,4:6]=dN-N/r
            end
            B[4,:]=B[1,:]+B[2,:]
            ke .+= weight*μ*(2B[1,:]*transpose(B[1,:])+2B[2,:]*transpose(B[2,:])+
                B[3,:]*transpose(B[3,:])-(2/3)*B[4,:]*transpose(B[4,:]))
            for c in 0:components-1
                ix=3c .+ (1:3); me[ix,ix] .+= weight*ρ*(N*transpose(N))
            end
            be .+= weight*B[4,:]*transpose(Q); qe .+= (weight/K)*(Q*transpose(Q))
        end
        for j in eachindex(uids),i in eachindex(uids)
            push!(I,uids[i]);push!(J,uids[j]);push!(mass,me[i,j]);push!(dev,ke[i,j])
            push!(coat,materials[e]==:coating ? ke[i,j] : 0.0)
        end
        for j in 1:2,i in eachindex(uids)
            push!(Ib,uids[i]);push!(Jb,qids[j]);push!(bv,be[i,j])
        end
        for j in 1:2,i in 1:2
            push!(Iq,qids[i]);push!(Jq,qids[j]);push!(qv,qe[i,j])
        end
    end
    (;M=sparse(I,J,mass,nu,nu),Kdev=sparse(I,J,dev,nu,nu),Kcoating=sparse(I,J,coat,nu,nu),
      B=sparse(Ib,Jb,bv,nu,np),Q=sparse(Iq,Jq,qv,np,np),nv)
end
