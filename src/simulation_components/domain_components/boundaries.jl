# Cardinal directions
export North, South, East, West
# Boundary of the model - types for each of the 4 walls (north, south, east, west)
export OpenBoundary, PeriodicBoundary, CollisionBoundary, MovingBoundary

"""
    North<:AbstractDirection

A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the northern boundary in a rectangular domain.
"""
struct North<:AbstractDirection end

#=
    _boundary_info_from_extent(North, FT, x0, xf, y0, yf)

Create a rectangular polygon representing the Northern boundary of a rectangular domain
defined by the x-extent (x0, xf) and the y-extent (y0, yf). The points of this boundary will
be of float type `FT`. If the length of the boundary x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`yf` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `yf` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _boundary_info_from_extent(::Type{North}, ::Type{FT}, x0, xf, y0, yf) where FT
    Δx, Δy = (xf - x0)/2, (yf - y0)/2
    poly = _make_bounding_box_polygon(FT, x0 - Δx, xf + Δx, yf, yf + Δy)
    return poly, FT(yf)
end

#= For floe interaction points (fpoints) that cross over the North boundary, zero-out all
forces in in direction not perpendicular to North boundary wall. =#
function _normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{North, FT}) where FT
    forces[fpoints[:, 2] .>= boundary.val, 1] .= FT(0.0)
    return
end

"""
    South<:AbstractDirection

A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the southern boundary in a rectangular domain.
"""
struct South<:AbstractDirection end

#=
    _boundary_info_from_extent(South, FT, x0, xf, y0, yf)

Create a rectangular polygon representing the Southern boundary of a rectangular domain
defined by the x-extent (x0, xf) and the y-extent (y0, yf). The points of this boundary will
be of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `y0` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `y0` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _boundary_info_from_extent(::Type{South}, ::Type{FT}, x0, xf, y0, yf) where FT
    Δx, Δy = (xf - x0)/2, (yf - y0)/2
    poly = _make_bounding_box_polygon(FT, x0 - Δx, xf + Δx, y0 - Δy, y0)
    return poly, FT(y0)
end

#= For floe interaction points (fpoints) that cross over the South boundary, zero-out all
forces in in direction not perpendicular to South boundary wall. =#
function _normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{South, FT}) where FT
    forces[fpoints[:, 2] .<= boundary.val, 1] .= FT(0.0)
    return
end

"""
    East<:AbstractDirection


A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the eastern boundary in a rectangular domain.
"""
struct East<:AbstractDirection end

#=
    _boundary_info_from_extent(East, FT, x0, xf, y0, yf)

Create a rectangular polygon representing the Eastern boundary of a rectangular domain
defined by the x-extent (x0, xf) and the y-extent (y0, yf). The points of this boundary will
be of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `xf` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `xf` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _boundary_info_from_extent(::Type{East}, ::Type{FT}, x0, xf, y0, yf) where FT
    Δx, Δy = (xf - x0)/2, (yf - y0)/2
    poly = _make_bounding_box_polygon(FT, xf, xf + Δx, y0 - Δy, yf + Δy)
    return poly, FT(xf)
end

#= For floe interaction points (fpoints) that cross over the East boundary, zero-out all
forces in in direction not perpendicular to East boundary wall. =#
function _normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{East, FT}) where FT
    forces[fpoints[:, 1] .>= boundary.val, 2] .= FT(0.0)
    return
end

"""
    West<:AbstractDirection


A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the western boundary in a rectangular domain.
"""
struct West<:AbstractDirection end

#=
    _boundary_info_from_extent(West, FT, x0, xf, y0, yf)

Create a rectangular polygon representing the Western boundary of a rectangular domain
defined by the x-extent (x0, xf) and the y-extent (y0, yf). The points of this boundary will
be of float type `FT`.  If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `x0` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `x0` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _boundary_info_from_extent(::Type{West}, ::Type{FT}, x0, xf, y0, yf) where FT
    Δx, Δy = (xf - x0)/2, (yf - y0)/2
    poly = _make_bounding_box_polygon(FT, x0 - Δx, x0, y0 - Δy, yf + Δy)
    return poly, FT(x0)
end

#= For floe interaction points (fpoints) that cross over the West boundary, zero-out all
forces in in direction not perpendicular to West boundary wall. =#
function _normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{West, FT}) where FT
    forces[fpoints[:, 1] .<= boundary.val, 2] .= FT(0.0)
    return
end

# Concrete subtype of AbstractBoundary - see documentation below
struct OpenBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    OpenBoundary{D, FT} <: AbstractBoundary{D, FT}

A concrete subtype of [`AbstractBoundary`](@ref) that removes a floe from the simulation if
any of the floe's vertices overlap with the `OpenBoundary`. This is meant to simulate the
floe floating out of the simulation. The floe is immeditaly removed as there might not be
ocean `u` and `v` values outside of the domian and thus we wouldn't know how to evolve the
floe's trajectory. This boundary is direction `D` ([`AbstractDirection`](@ref)) of the
domain and has field data of type `FT <: AbstractFloat`.

## _Fields_
- $BOUNDARY_POLY_DEF
- $VAL_DEF

Here is how to construct an `OpenBoundary`:

    OpenBoundary(D, [FT = Float64]; grid  = nothing, x0 = nothing, xf = nothing, y0 = nothing, yf = nothing)

The user must specify which [`AbstractDirection`](@ref) the boundary is. The user also has an opportunity to specify which
float type `Float64` or `Float32` will be used for field data. The user then must either provide
an grid object ([`AbstractRectilinearGrid`](@ref)) for the domain to align with the grid edges, or provide
four values (`x0`, `xf`, `y0`, and `yf`) that define the x and y-extents the user wants for the
domain.

!!! note
    If the user chooses to specify the domain extents, they still must be within the
    grid in order to make a valid [`Model`](@ref).

## _Positional arguments_
- $D_DEF
- $FT_DEF

## _Keyword arguments_
- $GRID_DEF
- $X0_DEF
- $XF_DEF
- $Y0_DEF
- $YF_DEF

!!! note
    The user must either provide a `grid` OR all four (4) of `x0`, `xf`, `y0`, and `yf`.

## _Examples_
- Defining a Northern `OpenBoundary` with Float64 (default type) data using the `grid` keyword.
```jldoctest
julia> g = RegRectilinearGrid(x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> OpenBoundary(North; grid = g)
OpenBoundary{North, Float64}
  ⊢polygon points are defined by the following set: (-250000.0, 750000.0), (750000.0, 500000.0), (750000.0, 750000.0), (-250000.0, 500000.0)
  ∟val is 500000.0
```
- Defining a Southern `OpenBoundary` with Float32 data using the `x0`, `xf`, `y0` and `yf` keywords.
```jldoctest
julia> OpenBoundary(South, Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5)
OpenBoundary{South, Float32}
  ⊢polygon points are defined by the following set: (750000.0f0, 0.0f0), (750000.0f0, -250000.0f0), (-250000.0f0, 0.0f0), (-250000.0f0, -250000.0f0)
  ∟val is 0.0
```
"""
function OpenBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing,
    x0 = nothing, xf = nothing, y0 = nothing, yf = nothing,
) where {D, FT}
    if !isnothing(grid)
        x0, xf, y0, yf = _get_grid_extent(grid)
    elseif isnothing(x0) || isnothing(xf) || isnothing(y0) || isnothing(yf)
        throw(ArgumentError("To create an OpenBoundary, either provide a grid or a x0, xf, y0, AND yf."))
    end
    poly, val = _boundary_info_from_extent(D, FT, x0, xf, y0, yf)
    return OpenBoundary{D, FT}(poly, val)
end

# Pretty printing for OpenBoundary showing key types and fields
function Base.show(io::IO, open_bound::OpenBoundary{D, FT}) where {D, FT}
    overall_summary = "OpenBoundary{$D, $FT}"
    points_summary, val_summary = show_boundary_poly_val_strings(open_bound)
    print(io, overall_summary, "\n",
        "  ⊢", points_summary, "\n",
        "  ∟", val_summary)
end

# Concrete subtype of AbstractBoundary - see documentation below
struct PeriodicBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    PeriodicBoundary <: AbstractBoundary

A concrete subtype of [`AbstractBoundary`](@ref) that moves a floe from one side of the
domain to the opposite side of the domain ([`North`](@ref) to [`South`](@ref) and
[`East`](@ref) to [`West`](@ref) and visa versa) when its centroid crosses the
`PeriodicBoundary`, bringing the floe back into the domain. Due to this behavior,
`PeriodicBoundary` pairs are required to form a valid domain.  This boundary is direction `D`
([`AbstractDirection`](@ref)) of the domain and has field data of type `FT <: AbstractFloat`.

## _Fields_
- $BOUNDARY_POLY_DEF
- $VAL_DEF

Here is how to construct an `PeriodicBoundary`:

    PeriodicBoundary(D, [FT = Float64]; grid  = nothing, x0 = nothing, xf = nothing, y0 = nothing, yf = nothing)

The user must specify which [`AbstractDirection`](@ref) the boundary is. The user also has an opportunity to specify which
float type `Float64` or `Float32` will be used for field data. The user then must either provide
an grid object ([`AbstractRectilinearGrid`](@ref)) for the domain to align with the grid edges, or provide
four values (`x0`, `xf`, `y0`, and `yf`) that define the x and y-extents the user wants for the
domain.

!!! note
    If the user chooses to specify the domain extents, they still must be within the
    grid in order to make a valid [`Model`](@ref).

## _Positional arguments_
- $D_DEF
- $FT_DEF

## _Keyword arguments_
- $GRID_DEF
- $X0_DEF
- $XF_DEF
- $Y0_DEF
- $YF_DEF

!!! note
    The user must either provide a `grid` OR all four (4) of `x0`, `xf`, `y0`, and `yf`.

## _Examples_
- Defining an Eastern `PeriodicBoundary` with Float32 data using the `grid` keyword.
```jldoctest
julia> g = RegRectilinearGrid(Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> PeriodicBoundary(East; grid = g)
PeriodicBoundary{East, Float64}
  ⊢polygon points are defined by the following set: (750000.0, -250000.0), (500000.0, -250000.0), (750000.0, 750000.0), (500000.0, 750000.0)
  ∟val is 500000.0
```
- Defining a Western `PeriodicBoundary` with Float64 data using the `x0`, `xf`, `y0` and `yf` keywords.
```jldoctest
julia> PeriodicBoundary(West, Float64; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5)
PeriodicBoundary{West, Float64}
  ⊢polygon points are defined by the following set: (0.0, -250000.0), (-250000.0, 750000.0), (0.0, 750000.0), (-250000.0, -250000.0)
  ∟val is 0.0
```
"""
function PeriodicBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing,
    x0 = nothing, xf = nothing, y0 = nothing, yf = nothing,
) where {D, FT}
    if !isnothing(grid)
        x0, xf, y0, yf = _get_grid_extent(grid)
    elseif isnothing(x0) || isnothing(xf) || isnothing(y0) || isnothing(yf)
        throw(ArgumentError("To create an PeriodicBoundary, either provide a grid or an x0, xf, y0, AND yf."))
    end
    poly, val = _boundary_info_from_extent(D, FT, x0, xf, y0, yf)
    return PeriodicBoundary{D, FT}(poly, val)
end

# Pretty printing for PeriodicBoundary showing key types and fields
function Base.show(io::IO, periodic_bound::PeriodicBoundary{D, FT}) where {D, FT}
    overall_summary = "PeriodicBoundary{$D, $FT}"
    points_summary, val_summary = show_boundary_poly_val_strings(periodic_bound)
    print(io, overall_summary, "\n",
        "  ⊢", points_summary, "\n",
        "  ∟", val_summary)
end

# Concrete type of AbstractBoundary - see documentation below
struct CollisionBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    CollisionBoundary <: AbstractBoundary

A concrete subtype of [`AbstractBoundary`](@ref) that calculates collision forces of a floe
against the boundary if any of the floe's vertices overlap with the `CollisionBoundary`.
This is meant to simulate any barrier that might stop a floe from flowing into or out of a
given region, like the edges of a cove. With this type of wall, a `CollisionBoundary` is
treated as an immovable, unbreakable floe for the purposes of calculations. This boundary is
direction `D` ([`AbstractDirection`](@ref)) of the domain and has field data of type
`FT <: AbstractFloat`.

## _Fields_
- $BOUNDARY_POLY_DEF
- $VAL_DEF

Here is how to construct an `CollisionBoundary`:

    CollisionBoundary(D, [FT = Float64]; grid  = nothing, x0 = nothing, xf = nothing, y0 = nothing, yf = nothing)

The user must specify which [`AbstractDirection`](@ref) the boundary is. The user also has an opportunity to specify which
float type `Float64` or `Float32` will be used for field data. The user then must either provide
an grid object ([`AbstractRectilinearGrid`](@ref)) for the domain to align with the grid edges, or provide
four values (`x0`, `xf`, `y0`, and `yf`) that define the x and y-extents the user wants for the
domain.

!!! note
    If the user chooses to specify the domain extents, they still must be within the
    grid in order to make a valid [`Model`](@ref).

## _Positional arguments_
- $D_DEF
- $FT_DEF

## _Keyword arguments_
- $GRID_DEF
- $X0_DEF
- $XF_DEF
- $Y0_DEF
- $YF_DEF

!!! note
    The user must either provide a `grid` OR all four (4) of `x0`, `xf`, `y0`, and `yf`.

## _Examples_
- Defining an Northern `CollisionBoundary` with Float64 data using the `grid` keyword.
```jldoctest
julia> g = RegRectilinearGrid(Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> CollisionBoundary(North; grid = g)
CollisionBoundary{North, Float64}
  ⊢polygon points are defined by the following set: (-250000.0, 750000.0), (750000.0, 500000.0), (750000.0, 750000.0), (-250000.0, 500000.0)
  ∟val is 500000.0
```
- Defining a Western `CollisionBoundary` with Float64 data using the `x0`, `xf`, `y0` and `yf` keywords.
```jldoctest
julia> CollisionBoundary(West, Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5)
CollisionBoundary{West, Float32}
  ⊢polygon points are defined by the following set: (0.0f0, -250000.0f0), (-250000.0f0, 750000.0f0), (0.0f0, 750000.0f0), (-250000.0f0, -250000.0f0)
  ∟val is 0.0
```
"""
function CollisionBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing,
    x0 = nothing, xf = nothing, y0 = nothing, yf = nothing,
) where {D, FT}
    if !isnothing(grid)
        x0, xf, y0, yf = _get_grid_extent(grid)
    elseif isnothing(x0) || isnothing(xf) || isnothing(y0) || isnothing(yf)
        throw(ArgumentError("To create an CollisionBoundary, either provide a grid or an x0, xf, y0, AND yf."))
    end
    poly, val = _boundary_info_from_extent(D, FT, x0, xf, y0, yf)
    return CollisionBoundary{D, FT}(poly, val)
end

# Pretty printing for CollisionBoundary showing key types and fields
function Base.show(io::IO, collision_bound::CollisionBoundary{D, FT}) where {D, FT}
    overall_summary = "CollisionBoundary{$D, $FT}"
    points_summary, val_summary = show_boundary_poly_val_strings(collision_bound)
    print(io, overall_summary, "\n",
        "  ⊢", points_summary, "\n",
        "  ∟", val_summary)
end

# Concrete type of AbstractBoundary - see documentation below
mutable struct MovingBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
    u::FT
    v::FT
end

"""
    MovingBoundary <: AbstractBoundary

A concrete subtype of [`AbstractBoundary`](@ref) that can provide a compressive or shear 
force on the floes within the domain by moving parallel or perpendicular to the domain.
The `MovingBoundary` calcuates the forces on the floe exactly the same as a
[`CollisionBoundary`](@ref), by acting as a immovable, unbreakable floe. This boundary is
direction `D` ([`AbstractDirection`](@ref)) of the domain and has field data of type
`FT <: AbstractFloat`.

If a `North` or `South` wall has a non-zero `v` velocity this will provide a compressive
stress as the floe moves towards the center of the domain along the vector `x = cx` where
`(cx, cy)` is the centroid of the domain (or decompresses with opposite sign velocities).
The `poly` and  `val` fields are updated to represent the movement of the domain at the
user-provided velocities, `u` and `v` m/s.

Alternatively, if a `North` or `South` wall has a non-zero `u` velocity this will provide a 
shear stress as the domain "moves" along the line `y = y0` (for `North`) or `y = yf`
(for South). In this case, this only changes the frictional forces and we do not need up
actually change the `poly` or `val` fields.

## _Fields_
- $BOUNDARY_POLY_DEF
- $VAL_DEF
- `u::FT`: boundary's u-velocity
- `v::FT`: boundary's v-velocity

Here is how to construct an `MovingBoundary`:

    MovingBoundary(D, [FT = Float64]; u = 0.0, v = 0.0, grid  = nothing, x0 = nothing, xf = nothing, y0 = nothing, yf = nothing)

The user must specify which [`AbstractDirection`](@ref) the boundary is. The user also has an opportunity to specify which
float type `Float64` or `Float32` will be used for field data. The user then must either provide
an grid object ([`AbstractRectilinearGrid`](@ref)) for the domain to align with the grid edges, or provide
four values (`x0`, `xf`, `y0`, and `yf`) that define the x and y-extents the user wants for the
domain.

!!! note
    If the user chooses to specify the domain extents, they still must be within the
    grid in order to make a valid [`Model`](@ref).

## _Positional arguments_
- $D_DEF
- $FT_DEF

## _Keyword arguments_
- `u::FT`: boundary's u-velocity
- `v::FT`: boundary's v-velocity
- $GRID_DEF
- $X0_DEF
- $XF_DEF
- $Y0_DEF
- $YF_DEF

!!! note
    The user must either provide a `grid` OR all four (4) of `x0`, `xf`, `y0`, and `yf`.
    If the user does not provide values for `u` or `v`, the boundary will not move.

## _Examples_
- Defining an Northern `MovingBoundary` with Float64 data using the `grid` keyword. Assigning u-velocity of 0.5m/s.
```jldoctest
julia> g = RegRectilinearGrid(Float32; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> MovingBoundary(North; u = 0.5, grid = g)
MovingBoundary{North, Float64}
  ⊢polygon points are defined by the following set: (-250000.0, 750000.0), (750000.0, 500000.0), (750000.0, 750000.0), (-250000.0, 500000.0)
  ⊢val is 500000.0
  ⊢u-velocity of 0.5 m/s
  ∟v-velocity of 0.0 m/s
```
- Defining a Southern `MovingBoundary` with Float32 data using the `x0`, `xf`, `y0` and `yf` keywords.
Assigning u-velocity of 0.3 m/s and v-velocity of 0.25 m/s
```jldoctest
julia> MovingBoundary(South, Float32; u = 0.3, v = 0.25, x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5)
MovingBoundary{South, Float32}
  ⊢polygon points are defined by the following set: (750000.0f0, 0.0f0), (750000.0f0, -250000.0f0), (-250000.0f0, 0.0f0), (-250000.0f0, -250000.0f0)
  ⊢val is 0.0
  ⊢u-velocity of 0.3 m/s
  ∟v-velocity of 0.25 m/s
```
"""
function MovingBoundary(
    ::Type{D},  # Directions
    ::Type{FT} = Float64;
    u = 0.0, v = 0.0, grid = nothing,
    x0 = nothing, xf = nothing, y0 = nothing, yf = nothing,
) where {D, FT}
    if !isnothing(grid)
        x0, xf, y0, yf = _get_grid_extent(grid)
    elseif isnothing(x0) || isnothing(xf) || isnothing(y0) || isnothing(yf)
        throw(ArgumentError("To create an MovingBoundary, either provide a grid or x0, xf, y0, AND yf."))
    end
    if u == v == 0
        @warn "MovingBoundary velocities are both zero. Boundary will not move. Use keywords u and v to set velocities."
    end
    poly, val = _boundary_info_from_extent(D, FT, x0, xf, y0, yf)
    return MovingBoundary{D, FT}(poly, val, u, v)
end

# Return MovingBoundary velocity as assigned by the user.
_get_velocity(boundary::MovingBoundary, _, _) = (boundary.u, boundary.v)

#= Move North/South MovingBoundary according to v-velocity and length of timestep by
changing the `poly` and `val` fields. =#
function _update_boundary!(boundary::MovingBoundary{D, FT}, Δt) where {D <: Union{North, South}, FT}
    (x0, xf), (y0, yf) = GI.extent(boundary.poly)
    Δy = boundary.v * Δt
    y0, yf = y0 + Δy, yf + Δy
    boundary.poly = _make_bounding_box_polygon(FT, x0, xf, y0, yf)
    boundary.val += Δy
    return
end

#= Move East/West MovingBoundary according to u-velocity and length of timestep by
changing the `poly` and `val` fields. =#
function _update_boundary!(boundary::MovingBoundary{D, FT}, Δt) where {D <: Union{East, West}, FT}
    (x0, xf), (y0, yf) = GI.extent(boundary.poly)
    Δx = boundary.u * Δt
    x0, xf = x0 + Δx, xf + Δx
    boundary.poly = _make_bounding_box_polygon(FT, x0, xf, y0, yf)
    boundary.val += Δx
    return
end

# Pretty printing for CollisionBoundary showing key types and fields
function Base.show(io::IO, moving_bound::MovingBoundary{D, FT}) where {D, FT}
    overall_summary = "MovingBoundary{$D, $FT}"
    points_summary, val_summary = show_boundary_poly_val_strings(moving_bound)
    u_velocity_summary = "u-velocity of $(moving_bound.u) m/s"
    v_velocity_summary = "v-velocity of $(moving_bound.v) m/s"
    print(io, overall_summary, "\n",
        "  ⊢", points_summary, "\n",
        "  ⊢", val_summary, "\n",
        "  ⊢", u_velocity_summary, "\n",
        "  ∟", v_velocity_summary)
end

# Union of all non-peridic boundary types to use as shorthand for dispatch.
const NonPeriodicBoundary = Union{OpenBoundary, CollisionBoundary, MovingBoundary}

const StationaryBoundary{D, FT} = Union{OpenBoundary{D, FT}, CollisionBoundary{D, FT}, PeriodicBoundary{D, FT}} where {D, FT}

# Return zero-velocities for StationaryBoundary.
_get_velocity(::StationaryBoundary{D, FT}, _, _) where {D, FT} =  (zero(FT), zero(FT))

# Do NOT update a StationaryBoundary during timestep.
_update_boundary!(::StationaryBoundary, Δt) = return

#= Boundaries across from one another (North and South OR East and West) are "periodically
compatible" if they are either both periodic, or if neither is periodic. This is because if
a floe exits a periodic boundary, it must be able to re-enter the opposite boundary to
fulfill its definition of periodic. This function is used to define valid domain walls.=#
_periodic_compat(::PeriodicBoundary, ::PeriodicBoundary) = true
_periodic_compat(::PeriodicBoundary, ::NonPeriodicBoundary) = false
_periodic_compat(::NonPeriodicBoundary, ::PeriodicBoundary) = false
_periodic_compat(::NonPeriodicBoundary, ::NonPeriodicBoundary) = true
