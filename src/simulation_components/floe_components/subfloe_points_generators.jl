export AbstractSubFloePointsGenerator, MonteCarloPointsGenerator, SubGridPointsGenerator

"""
    YourSubFloePointGenerator{FT} <: AbstractSubFloePointsGenerator{FT}

Abstract type for parameters determining generation of sub-floe points used for interpolation. The points generated using these parameters will be used to find
stresses on each floe from the ocean and the atmosphere within the coupling steps of the model. The two existing types are [`MonteCarloPointsGenerator`](@ref) and
[`SubGridPointsGenerator`](@ref). Both of these have different benefits. `MonteCarloPointsGenerator` guarentee that all floes have a somewhat similar number of subfloe points.
However, since these points are randomly placed, they are not neccesarily evenly spread out. There may not even be one point per model grid cell, which causes errors when
two-way coupling as stress from ice to ocean is calcualted with subfloe points. On the other hand with a `SubGridPointsGenerator`, each floe has a number of points
proportional to its area. However, these points are mainly evenly spaced and it is guarenteed that there is at least one per every model grid cell.
Therefore, you must use `SubGridPointsGenerator` when two-way coupling. 

## _API_

When a user creates a new subtype of AbstractSubFloePointsGenerator, the following functions must be created:

- `generate_subfloe_points(point_generator::MonteCarloPointsGenerator{FT}, poly, centroid, area, status, rng)`

This function that dispatches off of the subtype of `AbstractSubFloePointsGenerator` to generate the points for a given floe. It must
return the x-coordinate points, the y-coordiante points, and the floe's status. This should be the same as the input status, unless the
requested number of points cannot be generated, then the function should return `remove` (a type of `Status` tag). Points generated must all
be within a given floe, which is why the shape of the floe is passed as a polygon with the `poly` argument and why the floe's `centroid` and
area are also passed (they are used in existing implementations).
"""
abstract type AbstractSubFloePointsGenerator{FT<:AbstractFloat} end

# See documentation below
struct MonteCarloPointsGenerator{FT <: AbstractFloat} <: AbstractSubFloePointsGenerator{FT}
    npoints::Int
    ntries::Int
    err::FT

    function MonteCarloPointsGenerator{FT}(
        npoints,
        ntries,
        err,
    ) where {FT <: AbstractFloat}
        if npoints < 1
            throw(ArgumentError("Interpolation cannot be preformed with no \
            monte carlo points. Field npoints must be positive."))
        end
        if ntries < 1
            throw(ArgumentError("Monte carlo points cannot be generated without\
             trying at least once. Field ntried must be positive."))
        end
        if err < 0 || err > 1
            throw(ArgumentError("Field err must be between 0 and 1."))
        end
        return new{FT}(npoints, ntries, err)
    end
end

"""
    MonteCarloPointsGenerator{FT} <: AbstractSubFloePointsGenerator{FT}

Subtype of AbstractSubFloePointsGenerator that defines parameters needed to generate a set of random monte carlo points
within a given floe.

## Fields:

- `npoints::Int`: number of monte carlo points to generate within the floe's bounding box. Note that the floe will not end up with
    this many points as all points outside of the floe will be removed. (Default = 1000, probably much higher than needed!)
- `ntries::Int`: number of tries to generate a set of points within the floe that have a smaller error than `err` (Default = 100)
- `err::Real`: the percent of the floe that can not be covered by monte carlo points for it to be still be a valid set of subfloe points (Default = 0.1)

Here is how to construct an `MonteCarloPointsGenerator`:

    MonteCarloPointsGenerator([FT = Float64]; npoints = 1000, ntries = 100, err = 0.1)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- - `npoints::Int`: number of monte carlo points to generate within the floe's bounding box. 
- `ntries::Int`: number of tries to generate a set of points within the floe that have a smaller error than `err`
- `err::Real`: the percent of the floe that can not be covered by monte carlo points for it to be still be a valid set of subfloe points

## _Examples_

- Creating a default `MonteCarloPointsGenerator`
```jldoctest atmos
julia> MonteCarloPointsGenerator()
MonteCarloPointsGenerator{Float64}(1000, 10, 0.1)
```

- Creating a `MonteCarloPointsGenerator` of type Float32
```jldoctest atmos
julia> MonteCarloPointsGenerator(Float32; npoints = 100, ntries = 5, err = 0.05)
MonteCarloPointsGenerator{Float32}(100, 5, 0.05f0)
```
"""
MonteCarloPointsGenerator(::Type{FT} = Float64; npoints::Int = 1000, ntries::Int = 10, err = 0.1) where {FT <: AbstractFloat} =
    MonteCarloPointsGenerator{FT}(npoints, ntries, err)


"""
    generate_subfloe_points(
        point_generator
        poly,
        centroid,
        area,
        status,
        rng
    )

Generate monte carlo points centered on the origin within the floe according to
parameters defined in the point_generator argument.

## _Positional arguments_
    point_generator     <MonteCarloPointsGenerator> monte carlo point generator
    poly                <Polygon> Polygon representing floe shape
    centroid            <Matrix> floe's centroid
    area                <AbstractFloat> floe's area
    status              <Status> floe status (i.e. active, fuse in simulation)
    rng                 <AbstractRNG> random number generator to generate monte carlo points
## _Returns_
    x_sub_floe  <Vector{FT}> vector of sub-floe grid points x-coords within floe
    y_sub_floe  <Vector{FT}> vector of sub-floe grid points y-coords within floe
    status      <Status> floe's status post generation, changed to remove if generation is unsuccessful
"""
function generate_subfloe_points(
    point_generator::MonteCarloPointsGenerator{FT},
    poly,
    centroid,
    area,
    status,
    rng
) where {FT <: AbstractFloat}
    count = 1
    err = FT(1)
    mc_x = zeros(FT, point_generator.npoints)
    mc_y = zeros(FT, point_generator.npoints)
    mc_in = fill(false, point_generator.npoints)
    # Find bounding box
    poly = _translate_poly(FT, poly, -GI.x(centroid), -GI.y(centroid))::Polys{FT}
    (xmin, xmax), (ymin, ymax) = GI.extent(poly)
    Δx, Δy = xmax - xmin, ymax - ymin
    while err > point_generator.err
        if count > 10
            err = 0.0
            status.tag = remove
        else
            mc_x .= xmin .+ Δx * rand(rng, FT, point_generator.npoints)
            mc_y .= ymin .+ Δy * rand(rng, FT, point_generator.npoints)
            mc_in .= [coveredby_poly((mc_x[i], mc_y[i]), poly) for i in eachindex(mc_x)]
            err = abs(sum(mc_in)/point_generator.npoints * (Δx * Δy) - area)/area
            count += 1
        end
    end
    mc_x = mc_x[mc_in]
    mc_y = mc_y[mc_in]
    if isempty(mc_x)
        status.tag = remove
    end

    return mc_x, mc_y, status
end

# See documentation below
struct SubGridPointsGenerator{FT <: AbstractFloat} <: AbstractSubFloePointsGenerator{FT}
    Δg::FT

    function SubGridPointsGenerator{FT}(Δg) where {FT <: AbstractFloat}
        if Δg <= 0
            throw(ArgumentError("Field Δg must be positive as it is the width \
            and height value of the sub-floe grid cells."))
        end
        return new{FT}(Δg)
    end
end

"""
    SubGridPointsGenerator{FT} <: AbstractSubFloePointsGenerator{FT}


Subtype of [`AbstractSubFloePointsGenerator`](@ref) that defines parameters needed to
generate a set of points on a "subgrid" within the floe where the subgrid is a
regular rectilinar grid with cells of size `Δg` in both width and height.
The user can define how fine that grid should be in comparison with the model's grid.

## Fields:
- `Δg::FT`: regular rectilinar grid cell size for sub-floe points used for coupling/interpolation with ocean/atmosphere

!!! note
    If two-way coupling, `Δg` should be smaller than the grid's Δx and Δy so that there is at least one point in each
    grid cell that the floe occupies.

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `Δg::FT`: regular rectilinar grid cell size for sub-floe points used for coupling/interpolation with ocean/atmosphere (Default = Nothing)
- `grid::RegRectilinearGrid`: simulation's grid, which will be used to determine a reasonable `Δg`, along with `npoint_per_cell` if `Δg` isn't provided (Default = Nothing)
- `npoint_per_cell::Int`: number of points per sub-floe grid cell, where the grid is redefined to have width and height equal to the minimum of Δx and Δy

## _Returns_
    SubGridPointsGenerator with Δg defined to be the user defined value or the minimum of Δx and Δy over the number of desired points per grid cell

## _Examples_

- Creating a `SubGridPointsGenerator` with `Δg`
```jldoctest
julia> SubGridPointsGenerator(Δg = 1000)
SubGridPointsGenerator{Float64}(1000.0)
```

- Creating a `SubGridPointsGenerator` with `grid` and `npoint_per_cell`
```jldoctest atmos
julia> grid = RegRectilinearGrid(; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 50, Ny = 50);

julia> SubGridPointsGenerator(Float32; grid, npoint_per_cell = 10)
SubGridPointsGenerator{Float32}(707.1068f0)
```
"""
function SubGridPointsGenerator(::Type{FT} = Float64;
    Δg = nothing, grid = nothing, npoint_per_cell = nothing,
) where {FT <: AbstractFloat}
    Δg_val = if !isnothing(Δg)
        Δg
    elseif !isnothing(grid) & !isnothing(npoint_per_cell)
        min(grid.Δx, grid.Δy) / npoint_per_cell / sqrt(2)
    else
        throw(ArgumentError("To create a SubGridPointsGenerator, either Δg must be provided, or both a grid and npoint_per_cell must be provided."))
    end
    return SubGridPointsGenerator{FT}(Δg_val)
end

"""
    generate_subfloe_points(
        point_generator,
        poly,
        centroid,
        area,
        status,
        rng
    )

Generate evenly spaced points within given floe coordinates to be used for
coupling. If only one point falls within the floe, return the floe's centroid.
## _Positional arguments_
    point_generator     <MonteCarloPointsGenerator> monte carlo point generator
    poly                <Polygon> Polygon representing floe shape
    centroid            <Matrix> floe's centroid
    area                <AbstractFloat> floe's area
    status              <Status> floe status (i.e. active, fuse in simulation)
    rng                 <AbstractRNG> random number generator to generate monte
                            carlo points
## _Returns_
    x_sub_floe  <Vector{FT}> vector of sub-floe grid points x-coords within floe
    y_sub_floe  <Vector{FT}> vector of sub-floe grid points y-coords within floe
    status      <Status> floe's status post generation, changed to remove if generation is unsuccessful
"""
function generate_subfloe_points(
    point_generator::SubGridPointsGenerator{FT},
    poly,
    centroid,
    area,
    status,
    rng
) where {FT <: AbstractFloat}
    poly = _translate_poly(FT, poly, -GI.x(centroid), -GI.y(centroid))::Polys{FT}
    (xmin, xmax), (ymin, ymax) = GI.extent(poly)
    xpoints = Vector{FT}()
    ypoints = Vector{FT}()
    local x1, y1
    for (i, p) in enumerate(GI.getpoint(GI.getexterior(poly)))
        # Determine points on edges
        x2, y2 = GI.x(p), GI.y(p)
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        Δx, Δy = x2 - x1, y2 - y1
        l = sqrt((Δx)^2 + (Δy)^2)
        # Add current vertex
        push!(xpoints, x1)
        push!(ypoints, y1)
        # If distance between i and i+1 vertex is less than 2 sub-grid cells
        # but greater than one, add a point inbetween those two vertices
        x2_unshifted, y2_unshifted = x2, y2
        if l <= 2point_generator.Δg
            if l > point_generator.Δg
                push!(xpoints, x1 + Δx/2)
                push!(ypoints, y1 + Δy/2)
            end
        else  # The edge needs more points than the corners and midpoint
            if Δx == 0
                y1 +=  point_generator.Δg/2 * sign(Δy)
                y2 -= point_generator.Δg/2 * sign(Δy)
            elseif Δy == 0
                x1 += point_generator.Δg/2 * sign(Δx)
                x2 -= point_generator.Δg/2 * sign(Δx)
            else  # shift points to still be on the line
                m = Δy / Δx
                x_shift = sqrt(point_generator.Δg^2 / 4(1 + m^2))
                y_shift = m * x_shift
                x1 += x_shift
                x2 -= x_shift
                y1 += y_shift
                y2 -= y_shift
            end
            l = sqrt((x2 - x1)^2 + (y2 - y1)^2)
            n_edge_points = ceil(Int, l / point_generator.Δg) + 1
            append!(xpoints, range(x1, x2, length = n_edge_points))
            append!(ypoints, range(y1, y2, length = n_edge_points))
        end
        x1, y1 = x2_unshifted, y2_unshifted
    end
    # Add points in the interior of the floe
    n_xpoints = ceil(Int, (xmax - xmin) / point_generator.Δg)
    n_ypoints = ceil(Int, (ymax - ymin) / point_generator.Δg)
    x_interior_points = if n_xpoints < 3
        n_xpoints = 1
        FT(0):FT(0)  # polygon is centered at the origin
    else
        range(
            xmin + point_generator.Δg/2,
            xmax - point_generator.Δg/2,
            length = n_xpoints,
        )
    end
    y_interior_points = if n_ypoints < 3
        n_ypoints = 1
        FT(0):FT(0)
    else
        range(
            ymin + point_generator.Δg/2,
            ymax - point_generator.Δg/2,
            length = n_ypoints,
        )
    end
    x_sub_floe = repeat(x_interior_points, n_ypoints)
    y_sub_floe = repeat(y_interior_points, inner = n_xpoints)
    # Check which points are within the polygon and add to list
    in_floe = [coveredby_poly((x_sub_floe[i], y_sub_floe[i]), poly) for i in eachindex(x_sub_floe)]
    append!(xpoints, x_sub_floe[in_floe])
    append!(ypoints, y_sub_floe[in_floe])
    return xpoints, ypoints, status
end
