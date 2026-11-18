export Constants

const CONSTS_DEF = "`consts::Constants`: simulation's constants"

# See documentation below
@kwdef struct Constants{FT<:AbstractFloat}
    ρo::FT = 1027.0             # Ocean density
    ρa::FT = 1.2                # Air density
    Cd_io::FT = 3e-3            # Ice-ocean drag coefficent
    Cd_ia::FT = 1e-3            # Ice-atmosphere drag coefficent
    Cd_ao::FT = 1.25e-3         # Atmosphere-ocean momentum drag coefficient
    f::FT = 1.4e-4              # Ocean coriolis frequency
    turnθ::FT = 15π/180         # Ocean turn angle
    L::FT = 2.93e5              # Latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    ν::FT = 0.3                 # Poisson's ratio
    μ::FT = 0.2                 # Coefficent of friction
    E::FT = 6e6                 # Young's Modulus
end

"""
    Constants{FT}

Each simulation needs a list of physical constants to be used in various simulation calculations.
All constants have default values that will be used if the user does not provide an alternative.
    
## _Fields_

- `ρo::FT`: Ocean density (default = 1027.0)
- `ρa::FT`: Air density (default = 1.2)
- `Cd_io::FT`: Ice-ocean drag coefficent (default = 3e-3)
- `Cd_ia::FT`: Ice-atmosphere drag coefficent (default = 1e-3)
- `Cd_ao::FT`: Atmosphere-ocean momentum drag coefficient (default = 1.25e-3)
- `f::FT`: Ocean coriolis frequency (default = 1.4e-4)
- `turnθ::FT`: Ocean turn angle (default = 15π/180)
- `L::FT`: Latent heat of freezing [Joules/kg] (default = 2.93e5)
- `k::FT`: Thermal conductivity of surface ice[W/(m*K)] (default = 2.14)
- `ν::FT`: Poisson's ratio (default = 0.3)
- `μ::FT`: Coefficent of friction (default = 0.2)
- `E::FT`: Young's Modulus (default = 6e6)

Here is how to construct a `Constants`:

    Constants([FT = Float64]; kwargs...)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

!!! note
    Young's Modulus is usually calculated using the total floe area after floe initialization in the original Subzero code:
    `E = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))`

## _Examples_
- Creating `Constants` with a zero Coriolis frequency
```jldoctest
julia> Constants(; f = 0)
Constants{Float64}
  ⊢ ρo = 1027.0
  ⊢ ρa = 1.2
  ⊢ Cd_io = 0.003
  ⊢ Cd_ia = 0.001
  ⊢ Cd_ao = 0.00125
  ⊢ f = 0.0
  ⊢ turnθ = 0.2617993877991494
  ⊢ L = 293000.0
  ⊢ k = 2.14
  ⊢ ν = 0.3
  ⊢ μ = 0.2
  ⊢ E = 6.0e6
```

- Creating `Constants` with all default constants, but Float32 float-type
```jldoctest
julia> Constants(Float32)
Constants{Float32}
  ⊢ ρo = 1027.0
  ⊢ ρa = 1.2
  ⊢ Cd_io = 0.003
  ⊢ Cd_ia = 0.001
  ⊢ Cd_ao = 0.00125
  ⊢ f = 0.00014
  ⊢ turnθ = 0.2617994
  ⊢ L = 293000.0
  ⊢ k = 2.14
  ⊢ ν = 0.3
  ⊢ μ = 0.2
  ⊢ E = 6.0e6
```
"""
Constants(::Type{FT}, args...; kwargs...) where {FT <: AbstractFloat} =
    Constants{FT}(args...; kwargs...)

# Additional definition to allow to default FT = Float64
Constants(args...; kwargs...) = Constants{Float64}(args...; kwargs...)

# Pretty printing for Constants showing key dimensions
function Base.show(io::IO, consts::Constants{FT}) where FT
    overall_summary = "Constants{$FT}"
    consts_summary = ""
    for name in fieldnames(Constants)
        consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(consts, name)) * "\n")
    end
    print(io, overall_summary, "\n", consts_summary)
end
