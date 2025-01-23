"""
Types to hold parameters for simulation's physical process settings
"""

"""
    FloeSettings

Settings needed to create floes within the model.
- ρi is the density of all floes within the model
- min_floe_area is the minimum floe area within the model before removal
- min_floe_height is the minimum floe height within the model before removal
- max_floe_height is the maximum floe height within the model before the height
    can't increase any further
- min_aspect_ratio is the minimum ratio between the x-length and y-length of any
    floe prior to removal
- subfloe_point_generator is the method of subfloe point generation for each
    floe within the model
- stress_calculator is the method of calculating current stress of floe
"""
@kwdef struct FloeSettings{
    FT <: AbstractFloat,
    GT <: AbstractSubFloePointsGenerator{FT},
    CT <: AbstractStressCalculator{FT},
}
    ρi::FT = 920.0
    min_floe_area::FT = 1e6
    min_floe_height::FT = 0.1
    max_floe_height::FT = 10.0
    min_aspect_ratio::FT = 0.05
    maximum_ξ::FT = 1e-5
    subfloe_point_generator::GT = MonteCarloPointsGenerator()
    stress_calculator::CT = DecayAreaScaledCalculator()

    function FloeSettings{FT, GT, CT}(
        ρi,
        min_floe_area,
        min_floe_height,
        max_floe_height,
        min_aspect_ratio,
        maximum_ξ,
        subfloe_point_generator,
        stress_calculator,
    ) where {FT <: AbstractFloat, GT <: AbstractSubFloePointsGenerator{FT}, CT <: AbstractStressCalculator}
        if ρi < 0
            @warn "Ice density can't be negative. Resetting to default values of 920."
            ρi = FT(920)
        end
        if min_floe_area < 0
            @warn "Floe area can't be negative. Resetting minimum floe area to 0 m^2."
            min_floe_area = FT(0)
        end
        if min_floe_height < 0
            @warn "Floe height can't be negative. Resetting minimum floe area to 0."
            min_floe_height = FT(0)
        end
        if max_floe_height < 0
            @warn "Floe height can't be negative. Resetting to default of 10m."
            min_floe_height = FT(0)
        end
        if min_aspect_ratio < 0 || min_aspect_ratio > 1
            @warn "Aspect ratio must be between 0 and 1. Resetting to default of 0.05."
            min_aspect_ratio = FT(0.05)
        end
        if maximum_ξ < 0
            @warn "Maximum rotational velocity must be greater than 0. Resetting to default of 1e-5."
            min_aspect_ratio = FT(0.05)
        end
        new{FT, GT, CT}(
            ρi,
            min_floe_area,
            min_floe_height,
            max_floe_height,
            min_aspect_ratio,
            maximum_ξ,
            subfloe_point_generator,
            stress_calculator,
        )
    end

    FloeSettings(
        ρi,
        min_floe_area,
        min_floe_height,
        max_floe_height,
        min_aspect_ratio,
        maximum_ξ,
        subfloe_point_generator::GT,
        stress_calculator::CT,
    ) where {GT <: AbstractSubFloePointsGenerator, CT <: AbstractStressCalculator} = 
        FloeSettings{Float64, GT, CT}(
            ρi,
            min_floe_area,
            min_floe_height,
            max_floe_height,
            min_aspect_ratio,
            maximum_ξ,
            subfloe_point_generator,
            stress_calculator,
        )
end

"""
    FloeSettings(::Type{FT}; subfloe_point_generator, kwargs...)

A float type FT can be provided as the first argument of any FloeSettings
constructor. A FloeSettings of type FT will be created by passing all other
arguments to the correct constructor. 
"""
FloeSettings(
    ::Type{FT};
    subfloe_point_generator::GT = MonteCarloPointsGenerator(FT),
    stress_calculator::CT = DecayAreaScaledCalculator(FT),
    kwargs...,
) where {FT <: AbstractFloat, GT <: AbstractSubFloePointsGenerator, CT <: AbstractStressCalculator} =
    FloeSettings{FT, GT, CT}(;
        subfloe_point_generator,
        stress_calculator,
        kwargs...,
    )

"""
    CouplingSettings

Settings needed for coupling within the model.
If coupling_on is true, the model will be coupled with the simulation's ocean
and atmosphere. The Δt determines how many simulation timesteps between
calculating ocean and atmospheric forces on the floes. Δd number of buffer grid
cells on each side of floe for monte carlo interpolation and mc_n is the number
of monte carlo points to attempt to generage for each floe. If two_way_coupling_on is
true then the simulation calculates the stress the ice/atmosphere put on the
ocean. 
"""
@kwdef struct CouplingSettings
    coupling_on::Bool = true
    Δt::Int = 10
    Δd::Int = 1
    two_way_coupling_on::Bool = false

    function CouplingSettings(
        coupling_on,
        Δt,
        Δd,
        two_way_coupling_on,
    )
        if coupling_on && Δt < 0
            @warn "Coupling can't occur on a multiple of negative timesteps. \
                Turning coupling off."
            coupling_on = false
        end
        if !coupling_on && two_way_coupling_on
            @warn "Can't calculate stresses on ocean from ice and atmosphere \
                without coupling. Turning two_way_coupling_on off."
            two_way_coupling_on = false
        end
        if Δd < 0
            @warn "Can't complete interpolation of ocean and atmosphere forces \
                with a buffer of less than 0 grid cells. Setting Δd = 0."
            Δd = 0
        end
        new(
            coupling_on,
            Δt,
            Δd,
            two_way_coupling_on,
        )
    end
end

"""
    CollisionSettings{FT<:AbstractFloat}

Settings needed for collisions within the model. 
If collisions_on is true, collisions will occur, else they will not.
The floe_floe_max_overlap defines the percentage of overlap allowed between
floes before marking them for ridging/rafting.
The floe_domain_max_overlap defines the percentage of overlap allowed between
floes and the domain (collision boundaries and topography) before removing the
floe from the simulation. 
Both floe_floe_max_overlap and floe_domain_max_overlap should be between 0-1 and
if a value < 0 is given or a value > 1 is given when collisions_on is true they
will be set to 0 and 1 respectively.
"""
@kwdef struct CollisionSettings{FT<:AbstractFloat}
    collisions_on::Bool = true
    floe_floe_max_overlap::FT = 0.55
    floe_domain_max_overlap::FT = 0.75

    function CollisionSettings{FT}(
        collisions_on,
        floe_floe_max_overlap, 
        floe_domain_max_overlap,
    ) where {FT<:AbstractFloat}
        if collisions_on
            if floe_floe_max_overlap > 1
                @warn "The maximum collisin overlap between floes can't be \
                    greater than 1. Setting to 1."
                floe_floe_max_overlap = FT(1)
            elseif floe_floe_max_overlap < 0
                @warn "The maximum collisin overlap between floes can't be \
                    less than 0. Setting to 0."
                floe_floe_max_overlap = FT(0)
            end

            if floe_domain_max_overlap > 1
                @warn "The maximum collisin overlap between floes and the \
                    domain can't be greater than 1. Setting to 1."
                floe_domain_max_overlap = FT(1)
            elseif floe_domain_max_overlap < 0
                @warn "The maximum collisin overlap between floes and the \
                    domain can't be less than 0. Setting to 0."
                floe_domain_max_overlap = FT(0)
            end
        end
        new{FT}(
            collisions_on,
            floe_floe_max_overlap,
            floe_domain_max_overlap,
        )
    end
    CollisionSettings(
        collisions_on,
        floe_floe_max_overlap,
        floe_domain_max_overlap,
    ) = CollisionSettings{Float64}(
        collisions_on,
        floe_floe_max_overlap,
        floe_domain_max_overlap,
    )
end

"""
    CollisionSettings(::Type{FT}, kwargs...)

A float type FT can be provided as the first argument of any CollisionSettings
constructor. A CollisionSettings of type FT will be created by passing all
other arguments to the correct constructor. 
"""
CollisionSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    CollisionSettings{FT}(;kwargs...)

"""
    FractureSettings{CT<:AbstractFractureCriteria}

Settings needed for fractures within the model. 
If fractures_on is true, fractures will occur, else they will not.
The criteria defines which fracture criteria are used to determine which floes
to fracture. The Δt determines how many simulation timesteps between fracturing
floes. If deform_on is true, then the floe will be deformed around floe
primarily causing the fracture, identified by the largest overlap area on the
most recent set of collisions. Npieces denotes how many pieces to try to split a
fractured floe into - 3 is suggested value.
"""
@kwdef struct FractureSettings{CT<:AbstractFractureCriteria}
    fractures_on::Bool = false
    criteria::CT = NoFracture()
    Δt::Int = 0
    deform_on::Bool = false
    npieces::Int = 3

    function FractureSettings{CT}(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria}
        if fractures_on
            if Δt < 0
                @warn "Fracturing can't occur with negative timesteps. Turning \
                    fracturing off."
                fractures_on = false
            elseif criteria isa NoFracture
                @warn "Fracturing can't occur on with NoFracture criteria. \
                    Turning fracturing off."
                fractures_on = false
            elseif npieces < 2
                @warn "Fracturing can't occur on with npieces < 2 as this \
                    won't split floe. Turning fracturing off."
                fractures_on = false
            end
        end
        if !fractures_on && deform_on
            @warn "Deformation can't occur on without fracturing. Turning \
                deformation off."
            deform_on = false
        end
        new{CT}(fractures_on, criteria, Δt, deform_on, npieces)
    end
    FractureSettings(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria} = 
        FractureSettings{CT}(
            fractures_on,
            criteria,
            Δt,
            deform_on,
            npieces,
        )
end


"""
    struct SimplificationSettings{FT<:AbstractFloat}

If smooth_vertices_on is true then floe's with more vertices than
max_vertices will be simplified every Δt_smooth timesteps. The tolerance is the
Douglas–Peucker algorithm tolerance in (m).
"""
@kwdef struct SimplificationSettings{FT<:AbstractFloat}
    smooth_vertices_on::Bool = true
    max_vertices::Int = 30
    tol::FT = 100.0
    Δt_smooth::Int = 20

    function SimplificationSettings{FT}(
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    ) where {FT<:AbstractFloat}
        if smooth_vertices_on && Δt_smooth < 0
            @warn "Floe smoothing can't occur on a multiple of negative \
                timesteps. Turning floe simplification off."
            smooth_vertices_on = false
        end
        new{FT}(
            smooth_vertices_on,
            max_vertices,
            tol,
            Δt_smooth,
        )
    end
    
    SimplificationSettings(
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    ) = SimplificationSettings{Float64}(
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    )
end

"""
    SimplificationSettings(::Type{FT}, kwargs...)

A float type FT can be provided as the first argument of any
SimplificationSettings constructor. A SimplificationSettings of type FT will be
created by passing all other arguments to the correct constructor. 
"""
SimplificationSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    SimplificationSettings{FT}(;kwargs...)

"""
    RidgeRaftSettings{FT <: AbstractFloat}

Settings needed for ridging and rafting within the model. The fields have the
following meanings:
- ridge_raft_on: a boolean flag for if ridging and rafting should be turned on
    in the simulation
- Δt: multiple of timesteps during which ridging and rafting code will run
- ridge_probability: the probability a floe ridges with another floe/domain if
    it meets all other criteria
- raft_probability: the probability a floe rafts with another floe/domain if it
    meets all other criteria
- min_overlap_frac: the minimum overlap area fraction between a floe and another
    floe/domain for that floe to ridge or raft
- min_ridge_height: the minimum floe height to ridge with a floe/domain
- max_floe_ridge_height: the maximum floe height to ridge with another
    floe
- max_domain_rdige_height: maximum floe height to ridge with a domain element
- max_floe_raft_height: maximum floe height to raft with another floe
- max_domain_raft_height: maximum floe height to raft with a domain element
- domain_gain_probability: the probalility that a floe that rafts with a domain
    element keeps all of its mass (0) or if that mass is removed and lost to the
    domain element (1).
"""
@kwdef struct RidgeRaftSettings{FT<:AbstractFloat}
    ridge_raft_on::Bool = false
    Δt::Int = 0
    ridge_probability::FT = 0.95
    raft_probability::FT = 0.95
    min_overlap_frac::FT = 0.01
    min_ridge_height::FT = 0.2
    max_floe_ridge_height::FT = 5.0
    max_domain_ridge_height::FT = 1.25
    max_floe_raft_height::FT = 0.25
    max_domain_raft_height::FT = 0.25
    domain_gain_probability::FT = 1.0

    function RidgeRaftSettings{FT}(
        ridge_raft_on,
        Δt,
        ridge_probability,
        raft_probability,
        min_overlap_frac,
        min_ridge_height,
        max_floe_ridge_height,
        max_domain_ridge_height,
        max_floe_raft_height,
        max_domain_raft_height,
        domain_gain_probability,
    ) where {FT<:AbstractFloat}
        if ridge_raft_on && Δt < 0
            @warn "Ridging and rafting can't occur on a multiple of negative \
                timesteps. Turning ridging and rafting off."
            ridge_raft_on = false
        end
        if ridge_probability > 1
            @warn "Floes can't have a greater ridge probability than 1. \
            Setting ridge probability to 1."
            ridge_probability = FT(1)
        elseif ridge_probability < 0
            @warn "Floes can't have a smaller ridge probability than 1. \
            Setting ridge probability to 0."
            ridge_probability = FT(0)
        end
        if raft_probability > 1
            @warn "Floes can't have a greater raft probability than 1. \
            Setting ridge probability to 1."
            raft_probability = FT(1)
        elseif raft_probability < 0
            @warn "Floes can't have a smaller raft probability than 1. \
            Setting ridge probability to 0."
            raft_probability = FT(0)
        end
        if min_overlap_frac > 1
            @warn "Floes can't overlap more than 100%, so min_overlap_frac \
            can't exceed 1. Setting min_overlap_frac to 1."
            min_overlap_frac = FT(1)
        elseif min_overlap_frac < 0
            @warn "Floes can't overlap less than 0%, so min_overlap_frac \
            can't be less than 0. Setting min_overlap_frac to 0."
            min_overlap_frac = FT(0)
        end
        if domain_gain_probability > 1
            @warn "Floes can't have a greater domain_gain_probability than 1. \
            Setting domain_gain_probability to 1."
            domain_gain_probability = FT(1)
        elseif domain_gain_probability < 0
            @warn "Floes can't have a smaller domain_gain_probability than 0. \
            Setting domain_gain_probability to 0."
            domain_gain_probability = FT(0)
        end
        new{FT}(
            ridge_raft_on,
            Δt,
            ridge_probability,
            raft_probability,
            min_overlap_frac,
            min_ridge_height,
            max_floe_ridge_height,
            max_domain_ridge_height,
            max_floe_raft_height,
            max_domain_raft_height,
            domain_gain_probability,
        )
    end
    RidgeRaftSettings(
        ridge_raft_on,
        Δt,
        ridge_probability,
        raft_probability,
        min_overlap_frac,
        min_ridge_height,
        max_floe_ridge_height,
        max_domain_ridge_height,
        max_floe_raft_height,
        max_domain_raft_height,
        domain_gain_probability
    ) = RidgeRaftSettings{Float64}(
        ridge_raft_on,
        Δt,
        ridge_probability,
        raft_probability,
        min_overlap_frac,
        min_ridge_height,
        max_floe_ridge_height,
        max_domain_ridge_height,
        max_floe_raft_height,
        max_domain_raft_height,
        domain_gain_probability,
    )
end

"""
    RidgeRaftSettings(::Type{FT}, kwargs...)

A float type FT can be provided as the first argument of any RidgeRaftSettings
constructor. A RidgeRaftSettings of type FT will be created by passing all other
arguments to the correct constructor. 
"""
RidgeRaftSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    RidgeRaftSettings{FT}(;kwargs...)



"""
    WeldSettings{FT<:AbstractFloat}


Settings needed for welding within the model. The fields have the following
meanings:
- weld_on: a boolean flag for if welding should be turned on in the simulation
- Δts: a list of multiples of timesteps during which welding code will run,
    welding will be run at multiples of all elements, each with domain split
    into corresponding Nx and Ny values
- Nxs: a list of number of x-directional bins to split the domain into at
    corresponding timesteps
- Nys: a list of number of x-directional bins to split the domain into at
    corresponding timesteps
- min_weld_area: minimum area a weld can create for two floes to weld
- max_weld_area: maximum area a weld can create for two floes to weld
- welding_coeff: non-dimensional parameter, multiplied by ratio of overlap
    between two floes to original floe area to determin probability that a floe
    will merge. The larger this is, the more likely floes are to weld.
    Probability with 5% overlap is `welding_coeff * (0.05) > rand()`
"""
@kwdef struct WeldSettings{FT<:AbstractFloat}
    weld_on::Bool = false
    Δts::Vector{Int} = Vector{Int}()
    Nxs::Vector{Int} = Vector{Int}()
    Nys::Vector{Int} = Vector{Int}()
    min_weld_area::FT = 1e6
    max_weld_area::FT = 2e9
    welding_coeff::FT = 150
    function WeldSettings{FT}(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    ) where {FT<:AbstractFloat}
        if weld_on && (isempty(Δts) || any(Δts .≤ 0))
            @warn "Welding can't occur without any given timesteps or with \
            negative timesteps. Turning welding off."
            weld_on = false
        elseif any(Nxs .< 1) || any(Nys .< 1)
            @warn "Can't split the grid into less than one row or column. \
            Turning welding off." 
            weld_on = false
        elseif !(length(Δts) == length(Nxs) == length(Nys))
            @warn "Length of timestep multiple list (Δts) must match length of \
            grid split lists Nxs and Nys. Turning welding off." 
            weld_on = false
        end
        # Sort by largest to smallest timestep multiples
        order = reverse!(sortperm(Δts))
        Δts .= Δts[order]
        Nxs .= Nxs[order]
        Nys .= Nys[order]
        new{FT}(
            weld_on,
            Δts,
            Nxs,
            Nys,
            min_weld_area,
            max_weld_area,
            welding_coeff,
        )
    end
    WeldSettings(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    ) = WeldSettings{Float64}(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    )
end

"""
    WeldSettings(::Type{FT}, kwargs...)

A float type FT can be provided as the first argument of any WeldSettings
constructor. A WeldSettings of type FT will be created by passing all other
arguments to the correct constructor. 
"""
WeldSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    WeldSettings{FT}(; kwargs...)


# CORNERS::Bool = false           # If true, corners of floes can break
# PACKING::Bool = false           # If true, floe packing is enabled