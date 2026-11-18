export FloeSettings

const FLOE_SETTINGS_DEF = "`floe_settings::FloeSettings`: simulation's floe settings"

# See below of documentation
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

    FloeSettings{FT, GT, CT}

When you create a floe or a set of floes, you have the option to create a floe settings object. This set of settings controls certian floe fields and calculations.

## _Fields_
  - `ρi::FT`: floe's density (Default = 920.0 g/L)
  - `min_floe_area::FT`: minimum floe area (Default = 1e6 m^2)
  - `min_floe_height::FT`: minimum floe height (Default = 0.1 m)
  - `max_floe_height::FT`: maximum floe height (Default = 10.0 m)
  - `min_aspect_ratio::FT`: minimum ratio between floe x-length and y-length by maximum coordiante values (Default = 0.05)
  - `maximum_ξ::FT`: the absolute maximum rotational velocity a floe can reach before it is capped at maximum_ξ (Default = 1e-5 rad/s)
  - `subfloe_point_generator::GT`: subtype of [`AbstractSubFloePointsGenerator`](@ref), which generates floe's subfloe points.
        These points which determines the method of subfloe point generation is used for each floe (Default = MonteCarloPointsGenerator)
  - `stress_calculator::CT`: subtype of [`AbstractStressCalculator`](@ref), which generates the calculator for stress, which in turn determines
        the method of calculating current stress of floes during the simulation (Default = DecayAreaScaledCalculator)

If any of the minimum values are exceeded, a floe is removed in the course of the simulation. If any of the maximum values are reached,
the value is capped at the given value.

Here is how to construct a `FloeSettings`:

    FloeSettings([FT = Float64]; kwargs..)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `FloeSettings` 
```jldoctest
julia> FloeSettings()
FloeSettings{Float64, MonteCarloPointsGenerator{Float64}, DecayAreaScaledCalculator{Float64}}
  ⊢ ρi = 920.0
  ⊢ min_floe_area = 1.0e6
  ⊢ min_floe_height = 0.1
  ⊢ max_floe_height = 10.0
  ⊢ min_aspect_ratio = 0.05
  ⊢ maximum_ξ = 1.0e-5
  ⊢ subfloe_point_generator = MonteCarloPointsGenerator{Float64}(1000, 10, 0.1)
  ⊢ stress_calculator = DecayAreaScaledCalculator{Float64}(0.2, 0.0)
```

- Creating a Float32 `FloeSettings` with a zero minimum floe area and a `SubGridPointsGenerator` generator
```jldoctest
julia> subfloe_point_generator =  SubGridPointsGenerator(Float32; Δg = 1000);

julia> FloeSettings(Float32; min_floe_area = 0, subfloe_point_generator)
FloeSettings{Float32, SubGridPointsGenerator{Float32}, DecayAreaScaledCalculator{Float32}}
  ⊢ ρi = 920.0
  ⊢ min_floe_area = 0.0
  ⊢ min_floe_height = 0.1
  ⊢ max_floe_height = 10.0
  ⊢ min_aspect_ratio = 0.05
  ⊢ maximum_ξ = 1.0e-5
  ⊢ subfloe_point_generator = SubGridPointsGenerator{Float32}(1000.0f0)
  ⊢ stress_calculator = DecayAreaScaledCalculator{Float32}(0.2f0, 0.0f0)
```
"""
FloeSettings(
    ::Type{FT};
    subfloe_point_generator::GT = MonteCarloPointsGenerator(FT),
    stress_calculator::CT = DecayAreaScaledCalculator(FT),
    kwargs...,
) where {FT <: AbstractFloat, GT <: AbstractSubFloePointsGenerator, CT <: AbstractStressCalculator} =
    FloeSettings{FT, GT, CT}(;
        subfloe_point_generator = subfloe_point_generator,
        stress_calculator = stress_calculator,
        kwargs...,
    )

# Pretty printing for FloeSettings showing key dimensions
function Base.show(io::IO, floe_settings::FloeSettings{FT, GT, CT}) where {FT, GT, CT}
    overall_summary = "FloeSettings{$FT, $GT, $CT}"
    consts_summary = ""
    for name in fieldnames(FloeSettings)
        consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(floe_settings, name)) * "\n")
    end
    print(io, overall_summary, "\n", consts_summary)
end