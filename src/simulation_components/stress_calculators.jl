# # DecayAreaScaledCalculator and DamageStressCalculator

export AbstractStressCalculator, DecayAreaScaledCalculator, DamageStressCalculator

#= 
## What is a stress calculator? 

A stress calculator defines how the stress for a floe is calculated from its interactions
and how that stress is used to determine fracture. 

The `AbstractStressCalculator` is the abstract type, with `DecayAreaScaledCalculator`
currently implemented and the infrastructure for `DamageStressCalculator` setup for future
implementation.

Each of these concrete types represent different methods of calculating a floe's stress. A
`FloeSettings` struct requires a concrete subtype of `AbstractStressCalculator`. The
functions for calculating stress are then dispatched off of this calculator type. Those
functions currently are: `_update_stress_accum!` and `_scale_principal_stress!`. These
functions need to be extended for any new concrete subtype of `AbstractStressCalculator`.
=#

const STRESS_CALCULATOR_DEF = ": calculator to determine floe's stress calculations from floe interactions."

"""
    abstract type AbstractStressCalculator{FT <: AbstractFloat}

Abstract super type for stress calculators, which calculate a floe's stress from floe
interactions. The subtypes serve as dispatch types for the following three calculation
methods.

## API
The following methods must be implemented for all subtypes:
- `_update_stress_accum!(stress_calculator::AbstractStressCalculator{FT}, curr_stress::Matrix{FT} , floe::FloeType{FT})`
- `_scale_principal_stress!(stress_calculator::AbstractStressCalculator{FT}, σvals::Matrix{FT}, floe::FloeType{FT}, floe_settings::FloeSettings)`

`_update_stress_accum!` is called in the `calc_stress!` function and takes the stress
calculator, the `floe`'s instantatious stress at the current timestep, and the `floe` itself
and updates the `floe`'s `stress_accum` field, which is used when determining floe fracture
based off of stress. Within the function, other floe fields can be updated as needed.

`_scale_principal_stress!` is called within the `find_σpoint` function which is called
within the `determine_fractures` function. This function takes the stress calculator, the
`floe`'s accumulated stress in prinicpal stress space (`σvals = eigvals(stress_accum)`), the
`floe` itself, and the `floe_settings` and scales `σvals` by some values/ratio using
physical properties within `floe` and `floe_settings`. This is done to approximate changing
the polygon defining the floe's fracture criteria without having to redefine the fracture 
criteria for each floe. This is almost like a proxy for damage. 
"""
abstract type AbstractStressCalculator{FT <: AbstractFloat}  end

#= Default function to NOT scale stress before checking if stress point in principal stress
space is inside or outside of the fracture criteria polygon. Any calculators that do NOT
want to scale stress in this way can simply skip implementing this function and fall back
on this default. =#
_scale_principal_stress!(::AbstractStressCalculator, σvals, floe, floe_settings) = return

const DECAY_ARG_WARNING = "λ must be between 0 and 1. Resetting to default value of 0.2."
"""
    DecayAreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}

Type of AbstractStressCalculator that implements stress calculations by accumulating each
timestep of stress using a decay equation. The decay aspect increases importance placed on
new damage. The decay equation is as follows:

`stress_accum = stress_accum(1-λ) + stress_instant(λ)

The area-scaling part of the stress calculations comes into play when deciding if a floe
will be fractured. The accumulated stress can be scaled by a ratio of
(floe.area/min_floe_area)^α. By changing the α value, either larger or smaller floes
fracture more eassily. The scaled value will not be saved as this is equivalent to scaling
the simulation fracture criteria (morh's cone or hibler's ellipse, etc.), but it is less
computationally intensive. By leaving the default `α = 0`, this extra scaling will not
take place.

## Fields:
- `λ::AbstractFloat`: decay parameter used when calculating accumulated stress. Should be between 0 - 1.
- `α::AbstractFloat`: Adjusts ellipse in stress space by raising the ratio to the floe's area over the simulation minimum floe size to `α`.
## Note:
- `λ` is used in `_update_stress_accum!`, whereas α is used in `_scale_principal_stress!`.
"""
@kwdef struct DecayAreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    λ::FT = 0.2
    α::FT = 0.0

    function DecayAreaScaledCalculator{FT}(λ, α) where FT
        if λ < 0 || λ > 1
            @warn DECAY_ARG_WARNING
            λ = FT(0.2)
        end
        return new{FT}(λ, α)
    end
end

"""
    DecayAreaScaledCalculator(::Type{FT}; kwargs...)

A float type `FT` can be provided as the first argument of any
`DecayAreaScaledCalculator` constructor. A `DecayAreaScaledCalculator` of type `FT`
will be created by passing all other arguments to the correct constructor. 
"""
DecayAreaScaledCalculator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    DecayAreaScaledCalculator{FT}(args...; kwargs...)

"""
    DecayAreaScaledCalculator(; kwargs...)

If type isn't specified as the first argument, `DecayAreaScaledCalculator`(; kwargs...) will
be of type `Float64` and the correct constructor will be called with all other arguments.
"""
DecayAreaScaledCalculator(args...; kwargs...) = DecayAreaScaledCalculator{Float64}(args...; kwargs...)

#= Updating `floe`'s accumulated and instantanious stress using a decay equation given the
current stress at a timestep `curr_stress::Matrix{FT}`.=#
function _update_stress_accum!(stress_calculator::DecayAreaScaledCalculator, curr_stress, floe)
    λ = stress_calculator.λ
    floe.stress_accum .= (1 - λ) * floe.stress_accum + λ * curr_stress
    return
end

#= This function scales a floe's stress in pricipal stress space (eigenvalues of the
acucmulated stress) by `M = (floe_area/min_floe_area).^α ` as an alternative method of
adjusting the fracture criteria boundaries in principal stress space.=#
function _scale_principal_stress!(stress_calculator::DecayAreaScaledCalculator, σvals, floe, floe_settings)
    stress_calculator.α == 0 && return
    M = (floe.area/floe_settings.min_floe_area).^stress_calculator.α
    σvals .*= M
    return 
end

"""
    DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}

Type of AbstractStressCalculator that calculates stress with each timestep with
`damage * stress_instant`, as suggested by Mukund Gupta. This method could keep track of
damage directly within each floe using the `damage` floe field, perhapes as a value between
0-1, and rather than calculating an "accumulated stress" to store damage, as done in
`DecayAreaScaledCalculator`.

In this calculator, the floe's `damage` field keeps track of an explicit parameter.
## Fields:

## Note:
This method is not implemented and thus throws an error upon creation. The infrastructure
for the damage calculator is provided, but functions that depend on this calculator need to
be implemented as detailed in the `AbstractStressCalculator` API. 
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    function DamageStressCalculator{FT}() where FT
        throw(ErrorException("DamageStressCalculator isn't implemented yet!! Take a look at the issue/comments if you're interested in implementing."))
        return new{FT}()
    end
end

"""
    DamageStressCalculator(::Type{FT}; kwargs...)

A float type `FT` can be provided as the first argument of any
`DamageStressCalculator` constructor. A `DamageStressCalculator` of type `FT`
will be created by passing all other arguments to the correct constructor. 
"""
DamageStressCalculator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    DamageStressCalculator{FT}(args...; kwargs...)

"""
    DamageStressCalculator(; kwargs...)

If type isn't specified, `DamageStressCalculator`(; kwargs...) will be of type
`Float64` and the correct constructor will be called with all other arguments.
"""
DamageStressCalculator(args...; kwargs...) = DamageStressCalculator{Float64}(args...; kwargs...)

#= Updating `floe`'s damage parameter depending on the current floe's instantanious stress
and then update the floe's accumulated stress based on this damage. =#
function _update_stress_accum!(::DamageStressCalculator, curr_stress, floe)
    #= TODO: add in the changes to floe.damage and add any needed parameters to the
    DamageStressCalculator struct!! 
    floe.damage = ... =#
    floe.stress_accum .= floe.damage * curr_stress
    return
end