"""
Structs and functions used to define a Subzero model
"""

"""
    domain_in_grid(domain, grid)

Checks if given rectangular domain is within given grid and gives user a warning
if domain is not of maximum possible size given grid dimensions.
Inputs:
    domain      <RectangularDomain>
    grid        <AbstractRectilinearGrid>
Outputs:
    <Boolean> true if domain is within grid bounds, else false
"""
function domain_in_grid(domain::Domain, grid::AbstractRectilinearGrid)
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

"""
Model which holds grid, ocean, atmos structs, each with the same underlying
float type (either Float32 of Float64) and size. It also holds the domain
information, which includes the topography and the boundaries. It holds a 
StructArray of floe structs, again each relying on the same underlying float
type. Finally, it also holds the maximum floe id used thus far in the
simulation. This should be the length of the floes array at the beginning of the
run. 
"""
struct Model{
    FT<:AbstractFloat,             # float type
    GT<:AbstractRectilinearGrid{FT},          # grid type
    DT<:Domain{                    # domain type
        FT,
        <:AbstractBoundary,
        <:AbstractBoundary,
        <:AbstractBoundary,
        <:AbstractBoundary,
    },
    FLT<:StructArray{<:Floe{FT}},  # floe list type
}
    grid::GT
    ocean::Ocean{FT}
    atmos::Atmos{FT}
    domain::DT
    floes::FLT  # See floes.jl for floe creation

    function Model{FT, GT, DT, FLT}(
        grid::GT,
        ocean::Ocean{FT},
        atmos::Atmos{FT},
        domain::DT,
        floes::FLT,
    ) where {
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
        if !domain_in_grid(domain, grid)
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

    Model(
        grid::GT,
        ocean::Ocean{FT},
        atmos::Atmos{FT},
        domain::DT,
        floes::FLT,
    ) where {
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
    } = 
        Model{FT, GT, DT, FLT}(grid, ocean, atmos, domain, floes)
end