export TopographyElement  # Topographic element within domain
export initialize_topography_field  # Function to create topography field for user

# Concrete subtype of AbstractDomainElement - see documentation below
struct TopographyElement{FT}<:AbstractDomainElement{FT}
    poly::Polys{FT}
    centroid::Vector{FT}
    rmax::FT
end

"""
    TopographyElement{FT}<:AbstractDomainElement{FT}

A concrete subtype of [`AbstractDomainElement`](@ref) that represent topography in a given
simulation. These `TopographyElement` objects will be treated as be treated as islands or
coastline within the model in that they will not move or break due to floe interactions, but
they will affect floes that collide with them. 

Similar to [`Floe`](@ref) objects, we store the shape of each `TopographyElement` with a
`poly` field. We also store the `centroid` and the shape's maximum radius (`rmax`) in order
to calculate simple interaction checks with Floes.

## _Fields_
- $POLY_DEF
- $CENTROID_DEF
- $RMAX_DEF

Here is how to construct an `TopographyElement`:

    TopographyElement([FT = Float64]; poly::Polys)

The user must provide a `GeoInterface` polygon of the specific type defined by `Polys`.
The user also has an opportunity to specify which float type `Float64` or `Float32` will be
used for field data.

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- $POLY_DEF

## _Note_: The user should **NOT** be using this constructor. The user should create a topography
field using [`initialize_topography_field`](@ref) rather than creating individual topography
elements. This allows the abstraction of the use of the `StructArrays` package away from the
user.

## _Examples_
```jldoctest topography
julia> import GeoInterface as GI;

julia> poly = GI.Polygon([[(0.0, 0.0), (0.0, 1e3), (1e3, 1e3), (1e3, 0.0), (0.0, 0.0)]]);

julia> TopographyElement(Float64; poly)
TopographyElement{Float64}
  ⊢centroid is (500.0, 500.0) in meters
  ∟maximum radius is 707.1067811865476 meters
```

```jldoctest topography
julia> TopographyElement(Float32; poly)
TopographyElement{Float32}
  ⊢centroid is (500.0, 500.0) in meters
  ∟maximum radius is 707.1068 meters
```
"""
function TopographyElement(::Type{FT} = Float64; poly::Polys) where FT
    # Clean up polygon and calculate centroid and maximum radius
    poly = GO.ClosedRing()(poly)
    rmholes!(poly)
    centroid = collect(GO.centroid(poly, FT)) # TODO: Remove collect once type is changed
    rmax = calc_max_radius(poly, centroid, FT)
    return TopographyElement{FT}(poly, centroid, rmax)
end

# Return zero-velocities for topography elements 
_get_velocity(::TopographyElement{FT}, _, _) where {FT} =  (zero(FT), zero(FT))

# No forces should be zero-ed out in collidions with topography elements, unlike boundaries.
_normal_direction_correct!(_, _, ::TopographyElement) = return

# Pretty printing for TopographyElement showing key types and fields
function Base.show(io::IO, topo_element::TopographyElement{FT}) where FT
    overall_summary = "TopographyElement{$FT}"
    centroid_summary = "centroid is ($(GI.x(topo_element.centroid)), $(GI.y(topo_element.centroid))) in meters"
    rmax_summary = "maximum radius is $(topo_element.rmax) meters"
    print(io, overall_summary, "\n",
        "  ⊢", centroid_summary, "\n",
        "  ∟", rmax_summary)
end

"""
    initialize_topography_field([::Type{FT} = Float64]; polys = nothing, coords = nothing)

This function allows for easy initialization of a field of [`TopographyElement(@ref)s and
collects them into a `StructArray` so that they can be passed into a [`Model`](@ref).

This is the suggested way to create a topography field for a simulation. Do NOT construct
individual `TopographyElement` objects as that method does not correct the field as a whole
to ensure no topography polygons overlap and does not create a struct array that can be
passed to a `Model`.

The user can create a topography field by either passing a list of polygons or by passing
a list of coordiantes, which will then be made into polygons.

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- $POLY_LIST_DEF
- $COORDS_LIST_DEF

!!! note
    Topography field elements must not be intersecting, so the corrections within this
    function may lead to a different number of polygons than input to make a valid
    topography field.

## _Examples_
- Defining a topography field with coordinates
```jldoctest topo_field
julia> coords = [[[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.0)]], [[(10.0, 0.0), (10.0, 1.0), (11.0, 1.0), (10.0, 0.0)]]];

julia> initialize_topography_field(Float64; coords)
2-element TopographyField{Float64} list:
 TopographyElement{Float64}
  ⊢centroid is (0.3333333333333333, 0.6666666666666666) in meters
  ∟maximum radius is 0.74535599249993 meters
 TopographyElement{Float64}
  ⊢centroid is (10.333333333333334, 0.6666666666666666) in meters
  ∟maximum radius is 0.7453559924999301 meters
```

- Defining a topography field with polygons
```jldoctest topo_field
julia> import GeoInterface as GI;

julia> polys = [GI.Polygon(c) for c in coords];

julia> initialize_topography_field(Float32; polys)
2-element TopographyField{Float32} list:
 TopographyElement{Float32}
  ⊢centroid is (0.33333334, 0.6666667) in meters
  ∟maximum radius is 0.745356 meters
 TopographyElement{Float32}
  ⊢centroid is (10.333333, 0.6666667) in meters
  ∟maximum radius is 0.74535626 meters
```

- Creating an empty topography field without polys or coords
```jldoctest
julia> initialize_topography_field(Float64)
0-element TopographyField{Float64} list
```
"""
function initialize_topography_field(::Type{FT} = Float64; polys = nothing, coords = nothing) where FT
    # Make sure input given (if given) is turned into a list of polygons
    if isnothing(polys) && !isnothing(coords)
        polys = [make_polygon(c) for c in coords]
    elseif isnothing(polys)  # & isnothing(coords)
        return StructArray{TopographyElement{FT}}(undef, 0)
    end
    # Make sure list of polygons is non-intersecting --> these polygons will form topography
    topo_multipoly = GO.DiffIntersectingPolygons()(GI.MultiPolygon(polys))
    # Make each polygon from multipolygon into a topography element
    topo_field = StructArray{TopographyElement{FT}}(undef, GI.npolygon(topo_multipoly))
    for (i, poly) in enumerate(GI.getpolygon(topo_multipoly))
        topo_field[i] = TopographyElement(FT; poly)
    end
    return topo_field
end

"""
    TopographyField{FT} 

Alias for StructArray type with TopographyElement elements with data of type `FT`.
    
!!! note
    This is not designed to be used by users, but is useful for dispatch by developers. 
"""
const TopographyField{FT} = StructArray{TopographyElement{FT}} where FT

# Pretty printing for TopographyField showing key types and elements
function Base.showarg(io::IO, ::TopographyField{FT}, toplevel) where FT
    print(io, "TopographyField{$FT} list")
end
