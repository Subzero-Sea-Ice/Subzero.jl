export AbstractDomainElement  # All elements that make up the domain
export AbstractDirection  # # Cardinal directions for directional boundaries --> see boundaries.jl
export AbstractBoundary   # Used to create concrete boundary types --> see boundaries.jl

# Abstract type definitions and general fallback functions for all subtypes
"""
    YourDomainElement{FT} <: AbstractDomainElement{FT}

An abstract type for all of the element that create the shape and behavior of the [`Domain`](@ref):
the 4 boundary walls that make up the rectangular domain and the topography within the `Domain`. 

When building a [`Domain`](@ref), the user will create four (4) boundary elements that are
concrete types of [`AbstractBoundary`](@ref), which is itself a subtype of
`AbstractDomainElement`, and a list of [`TopographyElement`](@ref) objects, which are also a
concrete subtype of `AbstractDomainElement`.

## _API_
If the user wanted to add more domain elements, other than boundaries and topography, they
would create a new subtype of `AbstractDomainElement` and would need to write methods for
the following functions:
- `_get_velocity(element::AbstractDomainElement{FT}, x::FT, y::FT)`
- `_normal_direction_correct!(forces::Matrix{FT}, fpoints{Matrix{FT}}, element::AbstractDomainElement{FT})`

The `_get_velocity` function gets the velocity of a domain `element` at the point (`x`, `y`).
Given that currently implemented domain elements cannot rotate, the point at which the
velocity is measured is irrelevant. Furthermore, all domain elements but
[`MovingBoundary`](@ref) are stationary and thus have zero-velocity. This function is called
in the `calc_friction_forces` function.

The `_normal_direction_correct` function updates the `forces` provided to it and zero's out
any forces not in the normal direction, if applicable. Right now, this is used to ensure
that the elastic collision with boundaries only produces a normal foces on floes. For other 
domain elements (topography), this function does nothing.

!!! note
    - If the user is wants to create a new boundary type, they must also implement the \
    function in the [`AbstractBoundary`](@ref) API.
    - If the user wanted new elements in the [`Domain`](@ref), the user would need to change the \
    `Domain` struct.
"""
abstract type AbstractDomainElement{FT<:AbstractFloat} end

"""
    YourDirection <: AbstractDirection

Each domain within a Subzero.jl [`Model`](@ref) must have four (4) boundaries (subtypes
of [`AbstractBoundary`](@ref)) where each of these boundaries is parametrically typed by the
direction of the boundary. The user will first choose one of the four cardinal directions,
the subtypes of `AbstractDirection`:
- [`North`](@ref)
- [`South`](@ref)
- [`East`](@ref)
- [`West`](@ref)

This abstract type is not meant to be extended by the user, unless the user wants to move
away from a rectangular domain with assigned cardinal directions for each wall. This would
a more major redesign and the user should check out the [developer documentation]("devdocs.md").

## _API
- `_boundary_info_from_extent(D::Type{AbstractDirection}, FT::Type{AbstractFloat}, x0, xf, y0, yf)`

The `_boundary_info_from_extent` function finds the values for the `poly` and `val` fields of
a subtype of [`AbstractBoundary`](@ref) of direction `D` given a region of the extent defined
by `x0`, `xf`, `y0`, and `yf` given that these are the minimum and maximum `x` and `y` values
of the region desired. The returned extent will have values of type `FT`.
"""
abstract type AbstractDirection end

"""
    AbstractBoundary{D, FT} <: AbstractDomainElement{FT}

When running a Subzero simulation, each simulation has a rectangular domain that contains
the ice pack. This domain is made up for four boundaries, one in each of the cardinal directions:
[`North`](@ref), [`South`](@ref), [`East`](@ref), and [`West`](@ref). The user will pick
four boundary types, each concrete subtypes of `AbstractBoundary` and assign each of them a
direction, `D` and a float type `FT`, by typing them
[parametrically]("https://docs.julialang.org/en/v1/manual/types/#Parametric-Types"). 

The user can pick from the currently implemented concrete subtypes of `AbstractBoundary`
([`OpenBoundary`](@ref), [`CollisionBoundary`](@ref), [`PeriodicBoundary`](@ref), and
[`MovingBoundary`](@ref)) or implement their own.

For now, each boundary wall is assumed to be a rectangle. For correct behavior, the corners
of each of the boundary rectangles should overlap to make sure that no floes reach outside
of the domain without touching a boundary. These rectangles should be represented as
polygons with the `poly` field within each boundary. Furthermore, each boundary should have
a `val` field that represents either the line `y = val` (for `North` and `South` walls)
or `x = val` (for `East` and `West` walls) that denotes the line that marks the "inner-most"
edge of the domain. If an ice floe passes over the boundary's `val` line, the floe interacts
with the boundary.

```text
 ________________
|__|____val___|__| <- North coordinates include corners
|  |          |  |
|  |          |  | <- East and west coordinates ALSO include corners
|  |          |  |
```

For ease of use, each boundary should have a constructor option where the user can provide
an [`AbstractRectilinearGrid`](@ref) concrete type to define a boundary where `val` is at
the edge of the grid so that the boundaries form a border right around the grid.

## _API_
In addition to the below function, any new `AbstractBoundary` subtype must also implement
[`AbstractDomainElement`](@ref) API functions.

- `_update_boundary!(boundary::AbstractBoundary, Δt::Int)`
- `_periodic_compat(boundary1::AbstractBoundary, boundary2::AbstractBoundary)`

The `_update_boundary!` function updates a boundary's position (by changing the `val` and
`poly` fields) at every timestep. Currently, only `MovingBoundary` elements are updated, and
their update depends on the length of the simulation's timestep, `Δt`.

The `_periodic_compat` function determines if a pair of functions are "periodically compatible", 
that is, could they form a pair on opposite edges of a domain to create a set of functioning
periodic boundaries. It is called when creating a [`Domain`](@ref).
"""
abstract type AbstractBoundary{
    D<:AbstractDirection,
    FT,
}<:AbstractDomainElement{FT} end

# Helper function for boundary pretty printing for boundaries with `poly` and `val` fields
function show_boundary_poly_val_strings(boundary::AbstractBoundary)
    points = join(Set(GI.getpoint(boundary.poly))|>collect, ", ")
    points_summary = "polygon points are defined by the following set: $points"
    val_summary = "val is $(boundary.val)"
    return points_summary, val_summary
end

# Constants used in documentation for sub-type fields and keywords
const StaticQuadrilateral{FT} =  GI.Polygon{false,false, SA.SVector{1, GI.LinearRing{false, false, SA.SVector{5, Tuple{FT, FT}}, Nothing, Nothing}},Nothing,Nothing} where FT

const D_DEF = "`D::Type{<:AbstractDirection}`: subtype of AbstractDirection used to represent \\
if a boundary is the North, South, East, or West boundary of a simulation"

const BOUNDARY_POLY_DEF = "`poly::StaticQuadrilateral{FT}`: rectangular polygon representing \\
the shape and location of the boundary on a given timestep with points of float type `FT`."

const VAL_DEF = "`val::FT`: number of float type `FT` representing either the line `y = val` \\
(for `North` and `South` walls) or `x = val` (for `East` and `West` walls) that denotes the \\
line marking the inner-most edge of the domain."

const X0_DEF = "`x0::FT`: minimum x-value of domain extent (region floes will be contained within)"

const XF_DEF = "`xf::FT`: maximum x-value of domain extent (region floes will be contained within)"

const Y0_DEF = "`y0::FT`: minimum y-value of domain extent (region floes will be contained within)"

const YF_DEF = "`yf::FT`: maximum y-value of domain extent (region floes will be contained within)"
