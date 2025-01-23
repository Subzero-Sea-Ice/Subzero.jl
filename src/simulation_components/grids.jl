export AbstractRectilinearGrid, RegRectilinearGrid

# See CellFloes documentation below
struct CellFloes{FT<:AbstractFloat}
    floeidx::Vector{Int}
    Δx::Vector{FT}
    Δy::Vector{FT}
end

"""

    CellFloes{FT}

A `CellFloes` struct represents a single grid cell and accumulates the indices of
floes with area in that grid cell. This is used to accumulate the forces of floes in a grid
cell on the ocean below it. Due to the prevalence of periodic boundaries, the grid cell
represented by a `CellFloes` object are centered on grid points (translated by `Δx/2` and
`Δy/2` in the x and y directions for a [`RegRectilinearGrid`](@ref)), rather than on the
grid cells defined by the grid object itself. Floes are recorded with their index in the 
list of floes. Furthermore, in a model with periodic boundaries, some floes may be in
multiple grid cells on different edges of the domain if they pass through a periodic
boundary. In these cases, the floe "linked" to with its index must be translated by
a vector to get its "ghost" on the other side of the domain. This non-integer translation
data is of float type `FT`. `CellFloes` are used with [`CellStresses`](@ref) objects, which
aggregate the stress from each of the floes recorded in a given `CellFloes` object.

##  _Fields_
- `floeidx`: see keyword arguments
- `Δx`: see keyword arguments
- `Δy`: see keyword arguments

Here is how to construct a `CellFloes` object:

    CellFloes([FT = Float64]; floeidx = nothing, Δx = nothing, Δy = nothing)

## _Positional arguments_
- $FT_DEF
## _Keyword arguments_
- `floeidx::Vector{Int}`: vector of floe indicies in the list of model floes for floes with area in the grid
- `Δx::Vector{FT}`: vector of x-translations for a floe at the corresponding index of the
`floeidx` vector to be in the cell represented with the `CellFloes` object.
- `Δy::Vector{FT}`: vector of y-translations for a floe at the corresponding index of the
`floeidx` vector to be in the cell represented with the `CellFloes` object.

!!! note
    If no keyword arguments are provide by the user, an `CellFloes` object with empty
    fields will be created. This is the **standard useage** of these objects and they are added to
    during the coupling step. If keyword arguments are provided, then all three must be provided
    and each vector must be the same size.
"""
function CellFloes(::Type{FT} = Float64; floeidx = nothing, Δx = nothing, Δy = nothing) where {FT}
    if isnothing(floeidx) || isnothing(Δx) || isnothing(Δy)
        floeidx = Vector{Int}()
        Δx = Vector{FT}()
        Δy = Vector{FT}()
    else
        @assert length(floeidx) == length(Δx) == length(Δy) "A CellFloes object requires that all fields be vectors of the same length."
    end
    return CellFloes{FT}(floeidx, Δx, Δy)
end


# Empties each of the three vectors (`floeidx`, `Δx`, and `Δy`) within a CellFloes object.
function Base.empty!(cell::CellFloes)
    empty!(cell.floeidx)
    empty!(cell.Δx)
    empty!(cell.Δy)
    return
end

const GRID_DEF = "`grid::AbstractRectilinearGrid`: subtype of AbstractRectilinearGrid representing the grid used within a simulation"

"""
    YourGridType{FT} <: AbstractRectilinearGrid{FT}

Each simulation run with Subzero.jl must be run on a grid. A grid defines the points at
which the ocean and atmosphere hold velocity data, as well as various other tracers, at grid
points. The user must choose an existing subtype of `AbstractRectilinearGrid` or implement a new
subtype to create a simulation. 

Each grid implementation must define the number and dimensions of each grid cell. Right now,
this assumes that the ocean and the atmosphere are on the same grid and that the grid is
rectangular. We might not want this to be true in the future. Furthermore, each concrete
implementation of `AbstractRectilinearGrid` that will be used to two-way couple with an ocean
or an atmosphere must have a field called `floe_locations` that is a matrix of
[`CellFloes`](@ref), with one element for each grid point. The user should not worry about
`CellFloes`, this is up to the developer to make sure that their grid contains and populates
this field. Again, in the future, this might be related to the ocean and/or atmosphere
rather than the grid. For more information, or if you're interested in working on this, see
issue [#107](@ref).

## _API_
The following methods must be implemented for all subtypes:
- `_get_grid_extent(grid::AbstractRectilinearGrid)`: 

The `_get_grid_extent` function takes in a concrete subtype of `AbstractRectilinearGrid` and
returns the grid's minimum x-value, maximum x-value, minimum y-value, and maximum y-value.

Given the original code was written for [`RegRectilinearGrid`](@ref) objects, the dispatch
isn't fully implemented and other functions should be added to this list. If you were
interested in adding a new subtype of `AbstractRectilinearGrid`, see  issue [#108](@ref).
"""
abstract type AbstractRectilinearGrid{FT<:AbstractFloat} end

# Concrete subtype of AbstractRectilinearGrid - see below for documentation
struct RegRectilinearGrid{FT}<:AbstractRectilinearGrid{FT}
    Nx::Int
    Ny::Int
    Δx::FT
    Δy::FT
    x0::FT
    xf::FT
    y0::FT
    yf::FT
    floe_locations::Matrix{CellFloes{FT}}
end

"""
    RegRectilinearGrid{FT} <: AbstractRectilinearGrid{FT}

A concrete implementation of an [`AbstractRectilinearGrid`](@ref) that represents a tessellation of 2-dimensional
Euclidean space into `n`-by-`m` congruent rectangles. Fields that hold float data are of
type `FT`, a concrete subtype of `AbstractFloat`.

- `Nx::Int`: number of grid cells in the x-direction
- `Ny::Int`: number of grid cells in the y-direction
- `Δx::FT`: grid cell width
- `Δy::FT`: grid cell height
- `x0::FT`: value of first x grid line
- `xf::FT`: value of final x grid line
- `y0::FT`: value of first y grid line
- `yf::FT`: value of final y grid line
- `floe_locations::Matrix{CellFloes}`: `Nx + 1` by `Ny + 1` matrix of [`CellFloes`](@ref)

Here is how to construct a `RegRectilinearGrid`:

    RegRectilinearGrid([FT = Float64]; x0, xf, y0, yf, Nx = nothing, Ny = nothing, Δx = nothing, Δy = nothing)

## _Positional arguments_
- $FT_DEF
## _Keyword arguments_
- `x0::FT`: value of first x grid line
- `xf::FT`: value of final x grid line
- `y0::FT`: value of first y grid line
- `yf::FT`: value of final y grid line
- `Nx::Int`: number of grid cells in the x-direction
- `Ny::Int`: number of grid cells in the y-direction
- `Δx::FT`: grid cell width
- `Δy::FT`: grid cell height

!!! note
    The user MUST provide `x0`, `xf`, `y0`, and `yf`; the user then has the choice to
    provide `Nx` and `Ny` OR `Δx` and `Δy`. If provided `Δx` doesn't evenly divide length
    `xf-x0` or `Δy` doesn't evenly divide `yf-y0`, you won't get full size grid. The grid will
    be "trimmed" to the nearest full grid square in both directions.

## _Examples_

- Defining a `RegRectilinearGrid` using `Nx` and `Ny`.
```jldoctest
julia> grid = RegRectilinearGrid(; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 10)
RegRectilinearGrid{Float64}
  ⊢x extent (0.0 to 500000.0) with 20 grid cells of size 25000.0 m
  ∟y extent (0.0 to 500000.0) with 10 grid cells of size 50000.0 m
```
- Defining a `RegRectilinearGrid` using `Δx` and `Δy`.
```jldoctest
julia> grid = RegRectilinearGrid(Float32; x0 = -5e5, xf = 5e5, y0 = 0.0, yf = 5e5, Δx = 5e4, Δy = 5e4)
RegRectilinearGrid{Float32}
  ⊢x extent (-500000.0 to 500000.0) with 20 grid cells of size 50000.0 m
  ∟y extent (0.0 to 500000.0) with 10 grid cells of size 50000.0 m
```
- Error due to attemping to define a `RegRectilinearGrid` using `Δx` and `Ny`.
```jldoctest
julia> grid = RegRectilinearGrid(; x0 = -5e5, xf = 5e5, y0 = 0.0, yf = 5e5, Δx = 5e4, Ny = 10)
ERROR: ArgumentError: To create a RegRectilinearGrid, either provide Δx and Δy OR Nx and Ny.
[...]
```
"""
function RegRectilinearGrid(::Type{FT} = Float64;
    x0, xf, y0, yf,
    Nx = nothing, Ny = nothing,
    Δx = nothing, Δy = nothing,
) where FT <: AbstractFloat
    # ensure x and y endpoints are properly ordered
    if x0 > xf
        @warn "$x0 can't be larger than $xf. Swapping the order of the grid x-endpoints"
        xf, x0 = x0, xf
    end
    if y0 > yf
        @warn "$y0 can't be larger than $yf. Swapping the order of the grid y-endpoints"
        yf, y0 = y0, yf
    end
    # determine number and size of gridpoints based on provided arguments
    if !isnothing(Nx) && !isnothing(Ny)  # determine size of grid cells based on number
        Δx = (xf - x0) / Nx
        Δy = (yf - y0) / Ny
    elseif !isnothing(Δx) && !isnothing(Δy)  # determine number of grid cells based on size
        Nx = floor(Int, (xf - x0) / Δx)
        Ny = floor(Int, (yf - y0) / Δy)
        # note: round grid size if requested cell size doesn't go evenly into endpoints
        xf = x0 + Nx * Δx
        yf = y0 + Ny * Δy
    else
        throw(ArgumentError("To create a RegRectilinearGrid, either provide Δx and Δy OR Nx and Ny."))
    end
    # create cell floe list
    floe_locations = [CellFloes(FT) for _ in 1:(Nx + 1), _ in 1:(Ny + 1)]
    # create RegRectilinearGrid object
    return RegRectilinearGrid{FT}(Nx, Ny, Δx, Δy, x0, xf, y0, yf, floe_locations)
end

# Return the grid's minimum x-value, maximum x-value, minimum y-value, and maximum y-value
_get_grid_extent(grid::RegRectilinearGrid) = (grid.x0, grid.xf, grid.y0, grid.yf)

# Pretty printing for RegRectilinearGrid showing key dimensions
function Base.show(io::IO, grid::RegRectilinearGrid{FT}) where FT
    overall_summary = "RegRectilinearGrid{$FT}"
    x_summary = "x extent ($(grid.x0) to $(grid.xf)) with $(grid.Nx) grid cells of size $(grid.Δx) m"
    y_summary = "y extent ($(grid.y0) to $(grid.yf)) with $(grid.Ny) grid cells of size $(grid.Δy) m"
    print(io, overall_summary, "\n",
        "  ⊢", x_summary, "\n",
        "  ∟", y_summary)
end
