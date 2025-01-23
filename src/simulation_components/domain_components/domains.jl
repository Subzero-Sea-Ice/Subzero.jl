# Domain definition (combines 4 boundaries and topography)
export Domain

struct Domain{FT, NB, SB, EB, WB, TT}
    north::NB
    south::SB   
    east::EB
    west::WB
    topography::TT

    function Domain{FT, NB, SB, EB, WB, TT}(north::NB, south::SB, east::EB, west::WB, topography::TT) where {
        FT<:AbstractFloat,
        NB<:AbstractBoundary{North, FT},
        SB<:AbstractBoundary{South, FT},
        EB<:AbstractBoundary{East, FT},
        WB<:AbstractBoundary{West, FT},
        TT<:TopographyField{FT},
    }
        if !_periodic_compat(north, south)
            throw(ArgumentError("North and south boundary walls are not \
                periodically compatable as only one of them is periodic."))
        elseif !_periodic_compat(east, west)
            throw(ArgumentError("East and west boundary walls are not \
                periodically compatable as only one of them is periodic."))
        elseif north.val < south.val
            throw(ArgumentError("North boundary value is less than south \
                boundary value."))
        elseif east.val < west.val
            throw(ArgumentError("East boundary value is less than west \
                boundary value."))
        end
        new{FT, NB, SB, EB, WB, TT}(north, south, east, west, topography)
    end
end

"""
    Domain{FT, NB, SB, EB, WB, TT} 

A simulation `Domain` holds four (4) boundary elements (concrete subtypes of
[`AbstractBoundary`](@ref)) and a list (potentially empty) of [`TopographyElement`](@ref)
objects. There must be a wall typed by each cardinal direction, and each wall must hold data of the
same float type `FT`. Additionally, the [`North`](@ref) boundary should be "higher" (on the
y-axis) than the [`South`](@ref) boundary and the [`East`](@ref) should be futher to the
"right" (on the x-axis) than the [`West`](@ref) boundary. Aditionally, the boundary walls
should be overlapping as detailed in [`AbstractBoundary`](@ref).

Additionally, the set of walls must be periodically compatible. This means that pairs of
opposite boundaries (`North` and `South` AND `East` and `West`) both need to be periodic if
one of them is periodic. This is because if a floe exits a periodic boundary, it must be
able to re-enter the opposite boundary to fulfill its definition of periodic.

## _Fields_
- `north::NB`: Northern boundary where `NB <: AbstractBoundary{North, FT}`
- `south::SB`: Southern boundary where `SB <: AbstractBoundary{South, FT}`
- `east::EB`: Eastern boundary where `EB <: AbstractBoundary{East, FT}`
- `west::WB`: Western boundary where `WB <: AbstractBoundary{West, FT}`
- `topography::tT`: Field of topography elements where `TT <: TopographyField{FT}`

!!! note
    - All `FT` values above must be the same float type to form a valid domain.
    - The code depends on the boundaries forming a rectangle oriented along the cartesian grid. \
    Other shapes/orientations are not supported at this time. 

Here is how to construct an `MovingBoundary`:

    Domain(; north::NB, south::SB, east::EB, west::WB, topography::TT = nothing)

The user must provide the four (4) boundaries. They then have an option to provide topography.
If no topography is provided, an empty `TopographyField` will be created with the same `FT`
as the first of the boundaries (which should be shared amoung boundaries).

## _Keyword arguments_
- `north::NB`: domain's northern boundary
- `south::SB`: domain's southern boundary
- `east::EB`: domain's eastern boundary
- `west::WB`: domain's western boundary
- `topography::TT`: domain's topography field, if there is one, else an empty field is created

## _Examples_
- Creating a `Domain` with _NO_ `topography`
```jldoctest domain
julia> grid = RegRectilinearGrid(Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> north = OpenBoundary(North, Float64; grid);

julia> south = OpenBoundary(South, Float64; grid);

julia> east = PeriodicBoundary(East, Float64; grid);

julia> west = PeriodicBoundary(West, Float64; grid);

julia> Domain(; north, south, east, west)
Domain
  ⊢Northern boundary of type OpenBoundary{North, Float64}
  ⊢Southern boundary of type OpenBoundary{South, Float64}
  ⊢Eastern boundary of type PeriodicBoundary{East, Float64}
  ⊢Western boundary of type PeriodicBoundary{West, Float64}
  ∟0-element TopograpahyElement{Float64} list
```

- Creating a `Domain` with `topography`
```jldoctest domain
julia> import GeoInterface as GI;

julia> topo_polys = [GI.Polygon([[(1e4, 1e4), (1e4, 3e4), (3e4, 1e4), (1e4, 1e4)]]), GI.Polygon([[(8e4, 8e4), (8e4, 9e4), (9e4, 9e4), (9e4, 8e4), (8e4, 8e4)]])];

julia> topography = initialize_topography_field(; polys = topo_polys);

julia> Domain(; north, south, east, west, topography)
Domain
  ⊢Northern boundary of type OpenBoundary{North, Float64}
  ⊢Southern boundary of type OpenBoundary{South, Float64}
  ⊢Eastern boundary of type PeriodicBoundary{East, Float64}
  ⊢Western boundary of type PeriodicBoundary{West, Float64}
  ∟2-element TopograpahyElement{Float64} list
```
"""
function Domain(;
    north::NB, south::SB, east::EB, west::WB,  # need four boundary walls
    topography::TT = nothing,                  # optional topography input
) where {NB, SB, EB, WB, TT}
    FT = _get_boundary_float_type(north)
    if isnothing(topography)  # create an empty topography field
        return Domain(; north, south, east, west, topography = initialize_topography_field(FT))
    end
    return Domain{FT, NB, SB, EB, WB, TT}(north, south, east, west, topography)
end

# Helper function to access `FT` type from boundary
_get_boundary_float_type(::AbstractBoundary{<:AbstractDirection, FT}) where FT = FT

# Pretty printing for CollisionBoundary showing key types and fields
function Base.show(io::IO, domain::Domain{FT, NB, SB, EB, WB, TT}) where {FT, NB, SB, EB, WB, TT}
    overall_summary = "Domain"
    north_summary = "Northern boundary of type $NB"
    south_summary = "Southern boundary of type $SB"
    east_summary = "Eastern boundary of type $EB"
    west_summary = "Western boundary of type $WB"
    topo_summary = "$(length(domain.topography))-element TopograpahyElement{$FT} list"
    
    print(io, overall_summary, "\n",
        "  ⊢", north_summary, "\n",
        "  ⊢", south_summary, "\n",
        "  ⊢", east_summary, "\n",
        "  ⊢", west_summary, "\n",
        "  ∟", topo_summary)
end

#= =#
function get_domain_element(domain, idx)
    if idx == -1
        return domain.north
    elseif idx == -2
        return domain.south
    elseif idx == -3
        return domain.east
    elseif idx == -4
        return domain.west
    else
        topo_idx = -(idx + 4)
        return get_floe(domain.topography, topo_idx)
    end
end
