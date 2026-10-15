export RectangularBrick, SinusoidalBrick, specimen_outlines, specimen_area, straight_reference

# SI CAD features below 1 µm are outside this Gmsh/OCC adapter's contract.
const MIN_CAD_FEATURE = 1e-6

Base.@kwdef struct RectangularBrick
    length::Float64 = .017
    height::Float64 = .007
end

"""Correct STG: shifted indentation profiles, constant H-A over their common span. SI."""
Base.@kwdef struct SinusoidalBrick
    length::Float64 = .017
    height::Float64 = .007
    amplitude::Float64 = .0025
    periods::Int = 2
    end_margin::Float64 = .0004
    samples::Int = 80
end

function validate_shape(b::RectangularBrick)
    all(x -> isfinite(x) && x > 0,(b.length,b.height)) ||
        throw(ArgumentError("brick dimensions must be positive and finite"))
    min(b.length,b.height) >= MIN_CAD_FEATURE ||
        throw(ArgumentError("brick dimensions must resolve the 1 µm CAD minimum"))
    b
end

function validate_shape(b::SinusoidalBrick)
    all(x -> isfinite(x) && x > 0,(b.length,b.height)) ||
        throw(ArgumentError("brick dimensions must be positive and finite"))
    isfinite(b.amplitude) && 0 <= b.amplitude < b.height ||
        throw(ArgumentError("STG amplitude must be in [0,height)"))
    isfinite(b.end_margin) && 0 <= 2b.end_margin < b.length ||
        throw(ArgumentError("STG end margins leave no active span"))
    b.periods > 0 && b.samples >= 8 && b.samples % (2b.periods) == 0 ||
        throw(ArgumentError("STG samples must be >=8 and divisible by twice the period count"))
    (b.end_margin == 0 || b.end_margin >= MIN_CAD_FEATURE) &&
        b.height-b.amplitude >= MIN_CAD_FEATURE &&
        stg_active_span(b)/b.samples >= MIN_CAD_FEATURE ||
        throw(ArgumentError("STG edges/thickness must resolve the 1 µm CAD minimum; use exactly zero for no end margin"))
    b
end

stg_active_span(b::SinusoidalBrick) = (b.length-2b.end_margin)*2b.periods/(2b.periods+1)

function stg_bottom_depth(b::SinusoidalBrick,x)
    span=stg_active_span(b)
    b.end_margin <= x <= b.end_margin+span || return 0.0
    b.amplitude*(1-cos(2π*b.periods*(x-b.end_margin)/span))/2
end

function stg_top_depth(b::SinusoidalBrick,x)
    span=stg_active_span(b); shift=span/(2b.periods)
    start=b.end_margin; stop=start+span
    start+shift <= x <= stop+shift || return 0.0
    source_x=x <= stop ? x : start+(x-stop)
    clamp(b.amplitude-stg_bottom_depth(b,source_x),0.0,b.amplitude)
end

function brick_polygon(b::RectangularBrick)
    validate_shape(b)
    [(-b.length/2,-b.height/2),(b.length/2,-b.height/2),
     (b.length/2,b.height/2),(-b.length/2,b.height/2)]
end

function brick_polygon(b::SinusoidalBrick)
    validate_shape(b)
    span=stg_active_span(b); shift=span/(2b.periods)
    stop=b.end_margin+span
    bottom_x=merge_geometry_nodes(vcat(0.0,[b.end_margin+span*i/b.samples for i in 0:b.samples],b.length))
    top_x=merge_geometry_nodes(vcat(0.0,
        [b.end_margin+shift+span*i/b.samples for i in 0:b.samples],stop,b.length))
    bottom=[(x-b.length/2,stg_bottom_depth(b,x)-b.height/2) for x in bottom_x]
    top=[(x-b.length/2,b.height/2-stg_top_depth(b,x)) for x in reverse(top_x)]
    vcat(bottom,top)
end

function merge_geometry_nodes(values)
    sorted=sort!(values); merged=Float64[]
    tolerance=32eps(maximum(abs,sorted))
    for x in sorted
        (isempty(merged) || x-last(merged)>tolerance) && push!(merged,x)
    end
    merged
end

function local_outline(s)
    isnothing(s.shape) ? [(s.radius*cos(2π*i/64),s.radius*sin(2π*i/64)) for i in 0:63] :
        brick_polygon(s.shape)
end

"""Nominal outlines in SI. Shared by geometry checks and clients; meshing uses these brick vertices."""
function specimen_outlines(s)
    outline=local_outline(s)
    [[(p[1]+c[1],p[2]+c[2]) for p in outline] for c in s.centers]
end

function specimen_half_extents(s)
    isnothing(s.shape) ? (s.radius,s.radius) : (s.shape.length/2,s.shape.height/2)
end

specimen_min_thickness(s) = isnothing(s.shape) ? 2s.radius :
    s.shape isa SinusoidalBrick ? s.shape.height-s.shape.amplitude : s.shape.height

specimen_area(b::RectangularBrick) = (validate_shape(b); b.length*b.height)
specimen_area(b::SinusoidalBrick) = (validate_shape(b); b.length*b.height-b.amplitude*stg_active_span(b))
specimen_area(s) = length(s.centers)*(isnothing(s.shape) ? π*s.radius^2 : specimen_area(s.shape))

function polygon_signed_distance(p,vertices)
    inside=false; distance=Inf
    for i in eachindex(vertices)
        a,b=vertices[i],vertices[mod1(i+1,length(vertices))]
        dx,dy=b[1]-a[1],b[2]-a[2]
        t=clamp(((p[1]-a[1])*dx+(p[2]-a[2])*dy)/(dx^2+dy^2),0.0,1.0)
        distance=min(distance,hypot(p[1]-a[1]-t*dx,p[2]-a[2]-t*dy))
        if (a[2]>p[2]) != (b[2]>p[2]) &&
            p[1] < a[1]+(p[2]-a[2])*dx/dy
            inside=!inside
        end
    end
    inside ? -distance : distance
end

function distance_to_inclusion(s,p,c)
    isnothing(s.shape) && return hypot(p[1]-c[1],p[2]-c[2])-s.radius
    polygon_signed_distance((p[1]-c[1],p[2]-c[2]),brick_polygon(s.shape))
end

function straight_reference(s)
    s.shape isa SinusoidalBrick || throw(ArgumentError("straight reference requires an STG specimen"))
    b=s.shape
    remake_experiment(s;shape=RectangularBrick(length=b.length,height=specimen_area(b)/b.length))
end
