# Model definition (combines grid, domain, ocean, atmosphere, and floes)
export Model

const MODEL_DEF = "`model::Model`: simulation's model"

#= The model domain must fit within the model grid. Further warn users if the
domain is smaller than the grid as excess calculations will be carried out.
Ideally, grid and domain cover the same area. =#
function _domain_in_grid(domain::Domain, grid::AbstractRectilinearGrid)
    northval = domain.north.val
    southval = domain.south.val
    eastval = domain.east.val
    westval = domain.west.val
    if (northval <= grid.yf &&
        southval >= grid.y0 &&
        eastval <= grid.xf &&
        westval >= grid.x0)
        if (northval != grid.yf ||
            southval != grid.y0 ||
            eastval != grid.xf ||
            westval != grid.x0)
            @warn "At least one wall of domain is smaller than grid. This \
                could lead to unneeded computation. Consider making grid \
                smaller or domain larger."
        end 
        return true
    end
    return false
end

# See documentation below
struct Model{FT, GT, DT, FLT}
    grid::GT
    ocean::Ocean{FT}
    atmos::Atmos{FT}
    domain::DT
    floes::FLT

    function Model{FT, GT, DT, FLT}(grid::GT, ocean::Ocean{FT}, atmos::Atmos{FT}, domain::DT, floes::FLT) where {
        FT<:AbstractFloat,
        GT<:AbstractRectilinearGrid{FT},
        DT<:Domain{
            FT,
            <:AbstractBoundary,
            <:AbstractBoundary,
            <:AbstractBoundary,
            <:AbstractBoundary,
        },
        FLT<:StructArray{<:Floe{FT}},
    }
        if !_domain_in_grid(domain, grid)
            throw(ArgumentError("Domain does not fit within grid."))
        elseif (
            size(ocean.u) != size(atmos.u) ||
            size(ocean.v) != size(atmos.v) ||
            size(ocean.temp) != size(atmos.temp)
        )
            throw(ArgumentError("Ocean and atmosphere are not on the same grid.\
                This is not supported yet."))
        end
        if any(ocean.temp .< atmos.temp)
            @warn "In at least one grid cell the atmosphere temperature is \
            warmer than the ocean. This is not a situation in which the \
            thermodynamics are setup for right now."
        end
        new{FT, GT, DT, FLT}(grid, ocean, atmos, domain, floes)
    end
end

"""
    Model{FT, GT, DT, FLT} 

A simulation `Model` holds a concrete subtype of [`AbstractRectilinearGrid`](@ref),
a [`Domain`](@ref), an [`Ocean`](@ref), an [`Atmos`](@ref), and a list of [`Floe`](@ref) objects.
These are all of the physical pieces of the simulation. 

To sucessfully create a model, the domain must fit within the grid, the ocean and atmosphere must be the same dimensions
(this could be changed later once the code can handle different sized grids), and all of the model elements must be of the
same float type `FT`, either `Float64` or `Float32`.

## _Fields_
- `grid::GT`: model's grid where `GT <: AbstractRectilinearGrid{FT}`
- `ocean::Ocean{FT}`: model's ocean 
- `atmos::Atmos{FT}`: model's atmosphere
- `domain::DT`: model's domain where `DT <: Domain{FT, NB, SB, EB, WB, TT}`
- `floes::FLT `: models' list of floes where `FLT <: <:StructArray{<:Floe{FT}}`

!!! note
    - All `FT` values above must be the same float type to form a valid model.

Here is how to construct a `Model`:

    Model(; grid::GT, ocean::Ocean{FT}, atmos::Atmos{FT}, domain::DT, floes::FLT)

## _Keyword arguments_
- `grid::GT`: model's grid
- `ocean::Ocean{FT}`: model's ocean 
- `atmos::Atmos{FT}`: model's atmosphere
- `domain::DT`: model's domain
- `floes::FLT `: models' list of floes

## _Examples_
- Creating a `Model`
# ```jldoctest model
julia> using Random

julia> grid = RegRectilinearGrid(Float64; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> north = CollisionBoundary(North, Float64; grid);

julia> south = CollisionBoundary(South, Float64; grid);

julia> east = CollisionBoundary(East, Float64; grid);

julia> west = CollisionBoundary(West, Float64; grid);

julia> domain = Domain(; north, south, east, west);

julia> ocean = Ocean(Float64; u = 0.5, v = 0.25, temp = 0.0, grid);

julia> atmos = Atmos(Float64; u = 0.0, v = 0.1, temp = 0.0, grid);

julia> floe_settings = FloeSettings(Float64);

julia> floes = initialize_floe_field(Float64, 3, [0.5], domain, 0.25, 0; floe_settings, rng = Xoshiro(1));

julia> Model(; grid, domain, ocean, atmos, floes)
Model{Float64, ...}

 ⊢RegRectilinearGrid{Float64}
  ⊢x extent (0.0 to 500000.0) with 20 grid cells of size 25000.0 m
  ∟y extent (0.0 to 500000.0) with 20 grid cells of size 25000.0 m

 ⊢Domain
  ⊢Northern boundary of type CollisionBoundary{North, Float64}
  ⊢Southern boundary of type CollisionBoundary{South, Float64}
  ⊢Eastern boundary of type CollisionBoundary{East, Float64}
  ⊢Western boundary of type CollisionBoundary{West, Float64}
  ∟0-element TopograpahyElement{Float64} list

 ⊢Ocean{Float64}
  ⊢Vector fields of dimension (21, 21)
  ⊢Tracer fields of dimension (21, 21)
  ⊢Average u-velocity of: 0.5 m/s
  ⊢Average v-velocity of: 0.25 m/s
  ∟Average temperature of: 0.0 C

 ⊢Atmos{Float64}
  ⊢Vector fields of dimension (21, 21)
  ⊢Tracer fields of dimension (21, 21)
  ⊢Average u-velocity of: 0.0 m/s
  ⊢Average v-velocity of: 0.1 m/s
  ∟Average temperature of: 0.0 C

 ⊢Floe List:
  ⊢Number of floes: 3
  ⊢Total floe area: 1.375278018545777e11
  ∟Average floe height: 0.25
# ```
"""
Model(;
    grid::GT,
    ocean::Ocean{FT},
    atmos::Atmos{FT},
    domain::DT,
    floes::FLT,
) where {FT, GT, DT, FLT} = 
    Model{FT, GT, DT, FLT}(grid, ocean, atmos, domain, floes)

# Pretty printing for Model showing key dimensions
function Base.show(io::IO, model::Model{FT, GT, DT, FLT}; digits = 5) where {FT, GT, DT, FLT}
    overall_summary = "Model{$FT, ...}"
    print(io, overall_summary,
        "\n ⊢", model.grid,
        "\n ⊢", model.domain,
        "\n ⊢", model.ocean,
        "\n ⊢", model.atmos,
        "\n ∟", model.floes, 
    )
end