export AbstractFractureCriteria, NoFracture, HiblerYieldCurve, MohrsCone, FractureSettings

"""
    abstract type AbstractFractureCriteria

Abstract super type for fracture criteria, which determines if a floe fractures given its
physical condition. Right now, fracture criteria define a set of vertices that define a
shape in principal stress space. The minimum and maximum eigenvalues of a floe's stress field
are then a point within the principal stress space. If that stress point falls outside of
the criteria-verticies defined vertices, it is a stress big enough to fracture the floe.
Otherwise the floe will not be fractured. If another form of fracture criteria is desired,
some additional code may need to be pulled into helper fuctions within the fracture.jl file
and then dispatched off of using subtypes of AbstractFractureCriteria. This is encuarged if It
is helpful!

Right now, the subtypes serve as dispatch types for the following two methods.

## API
The following methods must be implemented for all subtypes:
- `__update_criteria!(criteria::AbstractFractureCriteria, floes::StructArray{<:Floe})`
- `_determine_fractures(criteria::AbstractFractureCriteria, floes::StructArray{<:Floe}, floe_settings::FloeSettings)`

`__update_criteria!` is called in the `fracture_floes!` function and takes the fracture
criteria and the current ice floe pack and updates the fracture criteria given the state
of the ice floe pack if needed.

`_determine_fractures` is called in the `fracture_floes!` and takes in the fracture
criteria, the current ice floe pack, and the floe settings and determines which of the ice
floes fracture. It returns a `Boolean` vector equal in length to the `floes` list. If the
`ith` element of the returned vector is `true` then the `i`th floe in `floes` will fracture.

Abstract type for fracture criteria. Each struct of this type must have a
vertices field representing the criteria in principal stress space. For a given
polygon, the minimum and maximum eigenvalues of its stress field will be its
location in principal stress space. If that stress point falls outside of
the criteria-verticies defined polygon it is a stress great enough to fracture
the floe. Otherwise the floe will not be fractured.
Each fracture criteria type must also have an _update_criteria! function defined
that is used to update the criteria each timestep. If the criteria does not need
to be updated, this function can be empty.

See the existing concrete subtypes: [`NoFracture`](@ref), [`MohrsCone`](@ref), and [`HiblerYieldCurve`](@ref).
"""
abstract type AbstractFractureCriteria end


#= Default function to NOT update the fracture criteria before determining floe fractures as
the simulation progresses. Any criteria that should not update as the simulation runs can
simply skip implementing this function and fall back on this default. =#
__update_criteria!(::AbstractFractureCriteria, _) = nothing

"""
    NoFracture<:AbstractFractureCriteria

Type of AbstractFractureCriteria representing when fracturing functionality is turned off.
If this is the type provided to the simulation's `FractureSettings`, then fractures will not
occur. As fracturing will not occur, this subtype depends on the default API functions.
"""
struct NoFracture<:AbstractFractureCriteria end

#= No floes are fractured when the fracture criteria is `NoFracture` and this function 
won't even be called due to given the `fractures_on` keyword is `false`. It is defined for
API completeness. =#
_determine_fractures(::NoFracture, _, _) = nothing

# Concrete subtype of AbstractFractureCriteria - see documentation below
mutable struct HiblerYieldCurve{FT<:AbstractFloat}<:AbstractFractureCriteria
    pstar::FT
    c::FT
    poly::Polys{FT}
end

#=
    _calculate_hibler(FT, floes, pstar, c)

Calculate Hibler's Elliptical Yield Curve as described in his 1979 paper
"A Dynamic Thermodynamic Sea Ice Model".

Uses the mean height of the current floe field and two tuning parameters to
determine the shape of the elliptical yield curve.

Hibler's paper says that: Both pstar and c relate the ice strength to the
ice thickness and compactness. c is determined so that 10% open water
reduces the strength substantially and pstar is considered a free parameter.
=#
function _calculate_hibler(::Type{FT}, mean_height, pstar, c) where FT
    compactness = 1  # Could be a user input with future development
    p = pstar*mean_height*exp(-c*(1-compactness))
    α_range = range(zero(FT), FT(2π), length = 100)
    a = p*sqrt(2)/2
    b = a/2
    ring_coords = [(a*cos(α), b*sin(α)) for α in α_range]
    ring_coords[end] = ring_coords[1] # make sure first and last element are exactly the same
    # TODO: eventually make with SVectors! 
    poly = GI.Polygon([ring_coords])
    return _move_poly(poly, -p/2, -p/2,  π/4)
end

"""
    HiblerYieldCurve{FT} <: AbstractFractureCriteria

Concrete subtype of AbstractFractureCriteria that calculates Hibler's Elliptical
Yield curve using parameters `pstar`, `c`, and the current floe field. This is an
elliptical yield curve that determines if a floe fractures based off if its stress
in principal stress space is inside or outside of that elliptical yield curve.

##  _Fields_
- `pstar::AbstractFloat`: parameter used to tune ellipse for optimal fracturing (Default = 2.25e5)
- `c::AbstractFloat`: parameter used to tune ellipse for optimal fracturing (Default = 20)
- `poly::Polys{FT}`: polygon that defines the yield curve in principal stress space

!!! note
    Based on Hibler's 1979 paper "A Dynamic Thermodynamic Sea Ice Model". Hibler's paper
    says that: Both `pstar` and `c` relate the ice strength to the ice thickness and compactness.
    `c` is determined so that 10% open water reduces the strength substantially and `pstar` is considered
    a free parameter.

Here is how to construct a `HiblerYieldCurve` object:

    HiblerYieldCurve([FT = Float64]; floes, pstar, c)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `floes::StructArray{Floes}`: models's list of floes
- `pstar::AbstractFloat`: parameter used to tune ellipse for optimal fracturing (Default = 2.25e5)
- `c::AbstractFloat`: parameter used to tune ellipse for optimal fracturing (Default = 20)

## _Examples_
- Creating default `HiblerYieldCurve` 
```jldoctest hibler_setup
julia> using Random

julia> grid = RegRectilinearGrid(Float64; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 20);

julia> north, south, east, west = CollisionBoundary(North, Float64; grid), CollisionBoundary(South, Float64; grid), CollisionBoundary(East, Float64; grid), CollisionBoundary(West, Float64; grid);

julia> domain = Domain(; north, south, east, west);

julia> floes = initialize_floe_field(Float64, 3, [0.5], domain, 0.25, 0; floe_settings = FloeSettings(Float64), rng = Xoshiro(1));

julia> hibler = HiblerYieldCurve(Float64; floes)
HiblerYieldCurve{Float64}
  ⊢ pstar: 225000.0
  ⊢ c: 20.0
  ⊢ yield curve area: 2.48338091663081e9
  ⊢ yield curve centroid: (-28125.0, -28125.0)
```

```jldoctest hibler_setup
julia> hibler = HiblerYieldCurve(Float32; floes, c = 24)
HiblerYieldCurve{Float32}
  ⊢ pstar: 225000.0
  ⊢ c: 24.0
  ⊢ yield curve area: 2.4833815e9
  ⊢ yield curve centroid: (-28124.994f0, -28124.998f0)
```
"""
function HiblerYieldCurve(::Type{FT} = Float64; floes, pstar = 2.25e5, c = 20.0) where FT
    vertices = _calculate_hibler(FT, mean(floes.height), pstar, c)
    return HiblerYieldCurve{FT}(pstar, c, vertices)
end

# syntactic sugar so previous versions of the code still work - not suggested! 
HiblerYieldCurve(floes, pstar = 2.25e5, c = 20.0) = HiblerYieldCurve(; floes, pstar, c)

#=
Update the Hibler Yield Curve vertices based on the current set of floes. The
criteria changes based off of the average height of the floes.
=#
function _update_criteria!(criteria::HiblerYieldCurve{FT}, floes) where FT
    criteria.poly = _calculate_hibler(
        FT,
        mean(floes.height),
        criteria.pstar,
        criteria.c
    )
    return
end

# Pretty printing for HiblerYieldCurve showing key dimensions
function Base.show(io::IO, curve::HiblerYieldCurve{FT}; digits = 5) where {FT}
    overall_summary = "HiblerYieldCurve{$FT}"
    pstar_summary = "  ⊢ pstar: " * string(curve.pstar)
    c_summary = "  ⊢ c: " * string(curve.c)
    area_summary = "  ⊢ yield curve area: " * string(round(area_poly(curve.poly, FT); digits))
    centroid_summary = "  ⊢ yield curve centroid: " * string(round.(centroid_poly(curve.poly, FT); digits))
    print(io, overall_summary, "\n", pstar_summary, "\n", c_summary, "\n", area_summary, "\n", centroid_summary)
end

# Concrete subtype of AbstractFractureCriteria - see documentation below
struct MohrsCone{FT<:AbstractFloat}<:AbstractFractureCriteria
    poly::Polys{FT}
end

#=
Creates polygon from vertex values for Mohr's Cone (triangle in 2D) in principal stress space
Equations taken from original version of Subzero written in MATLAB.
=#
function _calculate_mohrs(::Type{FT}, σ1, σ2, σ11, σ22) where FT
    # TODO: eventually make with SVectors! 
    points = [(σ1, σ2),  (σ11, σ22), (σ22, σ11), (σ1, σ2)]
    return GI.Polygon([points])
end

#=
Calculate Mohr's Cone coordinates in principal stress space and turn into a polygon
based off of coefficent of friction, uniaxial compressive strength, and the x-coordiante of
one vertex of the hibler cone (triangle in 2D)

Equations taken from original version of Subzero written in MATLAB.
=#
function _calculate_mohrs(
    ::Type{FT},
    q = 5.2,
    σc = 2.5e5,
    σ11 = -3.375e4,
) where FT
    σ1 = ((1/q) + 1) * σc / ((1/q) - q)
    σ2 = q * σ1 + σc
    σ22 = q * σ11 + σc
    return _calculate_mohrs(FT, -σ1, -σ2, -σ11, -σ22)
end

"""
    MohrsCone{FT} <: AbstractFractureCriteria

Concrete subtype of AbstractFractureCriteria that creates a conical yield curve in principal stress space.
This conical yield curve determines if a floe fractures based off if its stress in principal stress space
is inside or outside of that elliptical yield curve.

##  _Fields_
- `poly::Polys{FT}`: polygon that defines the yield curve in principal stress space

!!! note
    Based on concepts from Weiss, Jérôme, and Erland M. Schulson. "Coulombic faulting from the grain
    scale to the geophysical scale: lessons from ice." Journal of Physics D: Applied Physics 42.21 (2009): 214017.
    Equations taken from original version of Subzero written in MATLAB.

Here is how to construct a `MohrsCone` object:

    MohrsCone([FT = Float64]; kwargs...)


## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- `q::AbstractFloat`: based on the coefficient of internal friction (``µ_i``) by ``(μ_i^2 + 1)^(1/2) + μ_i^2`` (Default = 5.2)
- `σc::AbstractFloat`: uniaxial compressive strength (Default = 2.5e5)
- `σ11::AbstractFloat`: negative of the x-coordinate of one vertex of cone (triangle in 2D) and negative
    of the y-coordinate of adjacend vertex in principal stress space (Default = -3.375e4)
- `σ1::AbstractFloat`: x-coordiante of first point in cone (Default = nothing)
- `σ2::AbstractFloat`: y-coordiante of first point in cone (Default = nothing)
- `σ22::AbstractFloat`: y-coordinate of one vertex of cone and negative of the x-coordinate of
    adjacend vertex in principal stress space (Default = nothing)

!!! note
    The user must supply either `q`, `σc`, and `σ11` OR all four of `σ1`, `σ2`, `σ11`, and `σ22`. 

## _Examples_
- Creating default `MohrsCone` 
```jldoctest
julia> MohrsCone()
MohrsCone{Float64}
  ⊢ Points: (59523.80952, 59523.80952), (33750.0, -74500.0), (-74500.0, 33750.0), (59523.80952, 59523.80952)
```

- Creating Float32 `MohrsCone` 
```jldoctest
julia> MohrsCone(Float32)
MohrsCone{Float32}
  ⊢ Points: (59523.81f0, 59523.81f0), (33750.0f0, -74500.0f0), (-74500.0f0, 33750.0f0), (59523.81f0, 59523.81f0)
```
"""
function MohrsCone(::Type{FT} = Float64; q = 5.2, σc = 2.5e5, σ11 = -3.375e4,
    σ1 = nothing, σ2 = nothing, σ22 = nothing,
) where FT
    poly = if isnothing(σ1) || isnothing(σ2) || isnothing(σ22)
        _calculate_mohrs(FT, q, σc, σ11)
    else
        _calculate_mohrs(FT, -σ1, -σ2, -σ11, -σ22)
    end
    return MohrsCone{FT}(poly)
end

# syntactic sugar so previous versions of the code still work - not suggested! 
MohrsCone(args...) = MohrsCone(; args...)

# Mohr's cone is not time or floe dependent so it doesn't need to be updates.
function _update_criteria!(::MohrsCone, floes)
    return
end

# Pretty printing for MohrsCone showing key dimensions
function Base.show(io::IO, cone::MohrsCone{FT}; digits = 5) where {FT}
    overall_summary = "MohrsCone{$FT}"
    points_summary = "  ⊢ Points: "
    for (i, point) in enumerate(GI.getpoint(cone.poly))
        if i > 1
            points_summary *= ", "
        end
        points_summary *= string(round.(point; digits))
    end
    print(io, overall_summary, "\n", points_summary)
end

"""
    FractureSettings{CT}

Settings needed for fractures within the model. 

## _Fields_
- `fractures_on::Bool`: if true, floes will fracture, else they will not (Default = false)
- `criteria::CT`: defines which fracture criteria (subtype of [`AbstractFractureCriteria`](@ref)) are used to determine which floes to fracture
- `Δt::Int`: determines how many simulation timesteps between attempting floe fracture (Default = 0)
- `deform_on::Bool`: if true, then the floe will be deformed around floe primarily causing the fracture, identified by the largest overlap area on the most recent set of collisions (Default = False)
- `npieces::Int`: how many pieces to try to split a fractured floe into (Default = 3)

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `FractureSettings` 
```jldoctest
julia> FractureSettings()
FractureSettings{NoFracture}
  ⊢ fractures_on = false
  ⊢ criteria = NoFracture
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```

- Creating a Float32 `FractureSettings` with fractures on
```jldoctest
julia> mohrs_cone_criteria = MohrsCone(Float32);

julia> FractureSettings(; fractures_on = true, criteria = mohrs_cone_criteria)
FractureSettings{MohrsCone{Float32}}
  ⊢ fractures_on = true
  ⊢ criteria = MohrsCone{Float32}
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```
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

# Pretty printing for FractureSettings showing key dimensions
function Base.show(io::IO, fracture_settings::FractureSettings{CT}) where {CT}
    overall_summary = "FractureSettings{$CT}"
    consts_summary = ""
    for name in fieldnames(FractureSettings)
        if name == :criteria
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(typeof(getfield(fracture_settings, name))) * "\n")
        else
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(fracture_settings, name)) * "\n")
        end
    end
    print(io, overall_summary, "\n", consts_summary)
end