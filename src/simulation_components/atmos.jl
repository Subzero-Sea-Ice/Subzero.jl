export Atmos

# See documentation below
struct Atmos{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}

    function Atmos{FT}(u, v, temp) where FT
        if !(size(u) == size(v))
            throw(ArgumentError("One or more of the atmosphere vector fields \
             aren't the same dimension."))
        end
        new{FT}(u, v, temp)
    end
end

"""
    Atmos{FT}

Each simulation needs atmosphere/wind values on the grid to drive the ice floes. The `Atmos`
struct holds the needed fields to perform coupling of vector fields (u-velocity and 
v-velocity) and tracer fields (temperature).

Right now, all values are stored on grid points as defined by an [`AbstractRectilinearGrid`](@ref).
This should eventually change to a c-grid. If interested in this change, see issue [#30](@ref).

Coupling beyond velocity (for example thermodynamically) is not very well developed. If you
are interested in working on this see issue [#9](@ref). The `hflx_factor` field is related
to thermodynamic coupling.

## _Fields_
- `u::Matrix{FT}`: atmosphere/wind u-velocities in the x-direction on each grid point
- `v::Matrix{FT}`: atmosphere/wind v-velocities in the y-direction on each grid point
- `temp::Matrix{FT}`: atmosphere/wind temperature on each grid point

!!! note
    If a periodic boundary is used in the domain, the last grid cell in that direction will
    not be used as it is equivalent to the first grid cell. Thus, for all of these fields, the
    first and last value in the x and/or y direction should be equal if the east-west or
    north-south boundary pair are periodic respectively. In the future, it might be worth having
    the atmos re-sizing to not included this repeated point dispatching off of the boundary types.

!!! note
    For each of the fields, the rows represent the `x`-values of the grid, while the columns
    represent the `y`-values of the grid. This makes it easy to index into the grid for a
    point `(x, y)`, although it does mean that the fields provided don't "look" like the grid
    if plotted directly and instead need to be transposed. 

Here is how to construct an `Atmos`:

    Atmos([FT = Float64]; u, v, temp, grid = nothing)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- $U_DEF
- $V_DEF
- $TEMP_DEF
- $GRID_DEF

!!! note
    The user MUST provide `u`, `v`, and `temp`, although they have the option of these
    being `Real` values or `Matrix{<:Real}`. If the user chooses to provide `Real` values for any
    of those fields, then a `grid` IS required to create a constant field of the correct size.
    If user chooses to provide a field, they should be of size `(grid.Nx + 1, grid.Ny + 1)` and
    the user doesn't require a `grid` input.

## _Examples_
- Creating `Atmos` with constant `u`-velocity, `v-velocity`, and `temp`
```jldoctest atmos
julia> grid = RegRectilinearGrid(; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 10);

julia> Atmos(; u = 0.5, v = 0.25, temp = 0.0, grid)
Atmos{Float64}
  ⊢Vector fields of dimension (21, 11)
  ⊢Tracer fields of dimension (21, 11)
  ⊢Average u-velocity of: 0.5 m/s
  ⊢Average v-velocity of: 0.25 m/s
  ∟Average temperature of: 0.0 C
```
- Creating `Atmos` with user-defined `u` and `v` and a constant `temp` with `Float32` data
```jldoctest atmos
julia> import Random.Xoshiro;

julia> seeded_rng = Xoshiro(200);

julia> u_field = rand(seeded_rng, grid.Nx + 1, grid.Ny + 1);

julia> v_field = rand(seeded_rng, grid.Nx + 1, grid.Ny + 1);

julia> Atmos(Float32; u = u_field, v = v_field, temp = 0.0, grid)
Atmos{Float32}
  ⊢Vector fields of dimension (21, 11)
  ⊢Tracer fields of dimension (21, 11)
  ⊢Average u-velocity of: 0.4965465 m/s
  ⊢Average v-velocity of: 0.5340565 m/s
  ∟Average temperature of: 0.0 C
```
- Trying to create an `Atmos` with a constant field and NO grid
```jldoctest
julia> Atmos(; u = 0.2, v = 0.1, temp = 0.0)
ERROR: ArgumentError: To create a matrix from the constant value provided, you must provide a grid.
[...]
```
"""
function Atmos(::Type{FT} = Float64; u, v, temp, grid = nothing) where FT
    # Create field from provided values
    u_field = _get_val_field(FT, u, grid)
    v_field = _get_val_field(FT, v, grid)
    temp_field = _get_val_field(FT, temp, grid)
    return Atmos{FT}(u_field, v_field, temp_field)
end

# Pretty printing for Atmos showing key dimensions
function Base.show(io::IO, atmos::Atmos{FT}) where FT
    overall_summary = "Atmos{$FT}"
    vector_summary = "Vector fields of dimension $(size(atmos.u))"
    avg_atmos_u = "Average u-velocity of: $(mean(atmos.u)) m/s"
    avg_atmos_v = "Average v-velocity of: $(mean(atmos.v)) m/s"
    tracer_summary = "Tracer fields of dimension $(size(atmos.temp))"
    avg_atmos_temp = "Average temperature of: $(mean(atmos.temp)) C"
    print(io, overall_summary, "\n",
        "  ⊢", vector_summary, "\n",
        "  ⊢", tracer_summary, "\n",
        "  ⊢", avg_atmos_u, "\n",
        "  ⊢", avg_atmos_v, "\n",
        "  ∟", avg_atmos_temp)
end
