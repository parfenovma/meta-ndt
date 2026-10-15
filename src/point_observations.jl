export solid_displacement, solid_trace, harmonic_trace

"""Complex displacement `(ux, uy)` in metres at a solid point (SI coordinates).

Interpolate the existing P1 triangular solution, without meshing or solving.
Points outside the mesh, including traction-free voids, are rejected: absence
of solid is not a zero displacement measurement. No Gridap dependency is needed.
"""
function solid_displacement(r::SolidLensResult,point)
    ids,weights=solid_point_weights(r,point)
    (sum(w*r.displacement_x[i] for (i,w) in zip(ids,weights)),
     sum(w*r.displacement_y[i] for (i,w) in zip(ids,weights)))
end

function solid_point_weights(r,point)
    length(point)==2 && all(isfinite,point) || throw(ArgumentError("invalid solid point"))
    x,y=point
    for j in axes(r.triangles,2)
        a,b,c=r.triangles[:,j]
        ax,ay=r.nodes[:,a]; bx,by=r.nodes[:,b].-(ax,ay); cx,cy=r.nodes[:,c].-(ax,ay)
        det=bx*cy-by*cx
        iszero(det) && continue
        px,py=x-ax,y-ay
        v=(px*cy-py*cx)/det; w=(bx*py-by*px)/det; u=1-v-w
        min(u,v,w)>=-1e-9 || continue
        return ((a,b,c),(u,v,w))
    end
    throw(ArgumentError("point is not in the solid mesh (outside or in a void)"))
end

"""P1 displacement history `(times, ux, uy)` in seconds/metres at a solid point.

Requires `simulate_pulse(s;retain_history=true)`. Samples EVERY integration
step, never interpolates sparse animation frames or synthesizes a harmonic
signal. Querying a new point does not assemble or solve the FEM system.
"""
function solid_trace(r::SolidPulseResult,point)
    size(r.history_x,2)==length(r.times) ||
        throw(ArgumentError("pulse history not retained; use retain_history=true"))
    ids,weights=solid_point_weights(r,point)
    sample(h)=[sum(w*Float64(h[i,j]) for (i,w) in zip(ids,weights)) for j in eachindex(r.times)]
    (;times=r.times,ux=sample(r.history_x),uy=sample(r.history_y))
end

"""Samples of `Re(A * exp(-im*2π*f*t))`, with `f` in Hz and `t` in seconds.

This reconstructs a stationary harmonic response, not a transient propagation
or a new time-domain FEM solve. Output units are the units of the phasor `A`.
"""
function harmonic_trace(amplitude::Number,frequency::Real,times::AbstractVector)
    isfinite(amplitude) && isfinite(frequency) && frequency>0 && all(isfinite,times) ||
        throw(ArgumentError("harmonic phasor/times must be finite and frequency positive"))
    [real(amplitude)*cos(2π*frequency*t)+imag(amplitude)*sin(2π*frequency*t) for t in times]
end
