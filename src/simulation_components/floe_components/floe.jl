# Floe definition
export Floe, FILL_VALUE

const FLOE_DEF = "`floe::Floe`: singular floe within the simulation"

const FILL_VALUE = 1e32

# See documentation below
@kwdef mutable struct Floe{FT<:AbstractFloat}
    # Physical Properties -------------------------------------------------
    poly::Polys{FT}         # polygon that represents the floe's shape
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    rmax::FT                # distance of vertix farthest from centroid (m)
    moment::FT              # mass moment of intertia
    angles::Vector{FT}      # interior angles of floe in degrees
    # Monte Carlo Points ---------------------------------------------------
    x_subfloe_points::Vector{FT}        # x-coordinates for integration with ocean/atmos centered at origin
    y_subfloe_points::Vector{FT}        # y-coordinates for integration with ocean/atmos centered at origin
    # Velocity/Orientation -------------------------------------------------
    α::FT = 0.0             # floe rotation from starting position in radians
    u::FT = 0.0             # floe x-velocity
    v::FT = 0.0             # floe y-velocity
    ξ::FT = 0.0             # floe angular velocity
    # Status ---------------------------------------------------------------
    status::Status = Status() # floe is active in simulation
    id::Int = 0             # floe id - set to index in floe array at start of
                            #   sim - unique to all floes
    ghost_id::Int = 0       # ghost id - if floe is a ghost, ghost_id > 0
                            #   representing which ghost it is
                            #   if floe is not a ghost, ghost_id = 0
    parent_ids::Vector{Int} = Vector{Int}()  # if the floe was originally part
                            # of one or several floes, list parent ids
    ghosts::Vector{Int} = Vector{Int}()  # indices of ghost floes of given floe
    # Forces/Collisions ----------------------------------------------------
    fxOA::FT = 0.0          # force from ocean and atmos in x direction
    fyOA::FT = 0.0          # force from ocean and atmos in y direction
    trqOA::FT = 0.0         # torque from ocean and Atmos
    hflx_factor::FT = 0.0   # heat flux factor can be multiplied by floe height
                            #   to get the heatflux
    overarea::FT = 0.0      # total overlap with other floe
    collision_force::Matrix{FT} = zeros(1, 2)
    collision_trq::FT = 0.0
    interactions::Matrix{FT} = zeros(0, 7)
    num_inters::Int = 0
    stress_accum::Matrix{FT} = zeros(2, 2)
    stress_instant::Matrix{FT} = zeros(2, 2)
    strain::Matrix{FT} = zeros(2, 2)
    damage::FT = 0.0        # damage to the floe (to be used in new stress calculator)
    # Previous values for timestepping  -------------------------------------
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep y-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep angular-velocity
end

#=
Floe fields stored as plain arrays so that they can be moved to the GPU with `adapt`.
The first index of every array is the floe index. Ragged fields are padded:
- `poly[i, j, d]` is coordinate `d` of vertex `j` of floe `i`'s exterior ring, for
  `j in 1:n_points[i]` (the last vertex repeats the first); padding is `FILL_VALUE`.
- `interactions[i, k, c]` is column `c` of interaction `k` of floe `i`, for
  `k in 1:num_inters[i]`; padding is zero.
- `collision_force[i, d]`, `stress_accum[i, r, c]`, `stress_instant[i, r, c]`, and
  `strain[i, r, c]` hold the floe's small vector/matrix fields.
`flags` is a per-floe bit field that kernels can set to report events to the host.
=#
struct FixedWidthFloes{
    FT<:AbstractFloat,
    AT3<:AbstractArray{FT, 3},
    AT2<:AbstractArray{FT, 2},
    AT1<:AbstractArray{FT, 1},
    IT1<:AbstractArray{Int32, 1},
    UT1<:AbstractArray{UInt8, 1},
}
    # Shape
    n_points::IT1
    poly::AT3
    centroid::AT2
    area::AT1
    height::AT1
    mass::AT1
    moment::AT1
    # Velocity/orientation
    α::AT1
    u::AT1
    v::AT1
    ξ::AT1
    # Forces/collisions
    fxOA::AT1
    fyOA::AT1
    trqOA::AT1
    hflx_factor::AT1
    collision_force::AT2
    collision_trq::AT1
    num_inters::IT1
    interactions::AT3
    stress_accum::AT3
    stress_instant::AT3
    strain::AT3
    damage::AT1
    # Previous values for timestepping
    p_dxdt::AT1
    p_dydt::AT1
    p_dudt::AT1
    p_dvdt::AT1
    p_dξdt::AT1
    p_dαdt::AT1
    # Events reported by kernels
    flags::UT1
end

Adapt.@adapt_structure FixedWidthFloes

# Copy of floes as a FixedWidthFloes on the CPU. Does not share memory with floes.
function FixedWidthFloes(floes::StructArray{<:Floe{FT}}) where FT
    nfloes = length(floes)
    n_points = Int32[GI.npoint(GI.getexterior(p)) for p in floes.poly]
    max_points = isempty(floes) ? 0 : Int(maximum(n_points))
    poly = fill(FT(FILL_VALUE), nfloes, max_points, 2)
    centroid = Matrix{FT}(undef, nfloes, 2)
    collision_force = Matrix{FT}(undef, nfloes, 2)
    num_inters = Int32.(floes.num_inters)
    max_inters = isempty(floes) ? 0 : Int(maximum(num_inters))
    interactions = zeros(FT, nfloes, max_inters, 7)
    stress_accum = Array{FT}(undef, nfloes, 2, 2)
    stress_instant = Array{FT}(undef, nfloes, 2, 2)
    strain = Array{FT}(undef, nfloes, 2, 2)
    for i in 1:nfloes
        for (j, point) in enumerate(GI.getpoint(GI.getexterior(floes.poly[i])))
            poly[i, j, 1], poly[i, j, 2] = get_tuple_point(point, FT)
        end
        for d in 1:2
            centroid[i, d] = floes.centroid[i][d]
            collision_force[i, d] = floes.collision_force[i][d]
        end
        interactions[i, 1:num_inters[i], :] .= @view floes.interactions[i][1:num_inters[i], :]
        stress_accum[i, :, :] .= floes.stress_accum[i]
        stress_instant[i, :, :] .= floes.stress_instant[i]
        strain[i, :, :] .= floes.strain[i]
    end
    return FixedWidthFloes(
        n_points,
        poly,
        centroid,
        copy(floes.area),
        copy(floes.height),
        copy(floes.mass),
        copy(floes.moment),
        copy(floes.α),
        copy(floes.u),
        copy(floes.v),
        copy(floes.ξ),
        copy(floes.fxOA),
        copy(floes.fyOA),
        copy(floes.trqOA),
        copy(floes.hflx_factor),
        collision_force,
        copy(floes.collision_trq),
        num_inters,
        interactions,
        stress_accum,
        stress_instant,
        strain,
        copy(floes.damage),
        copy(floes.p_dxdt),
        copy(floes.p_dydt),
        copy(floes.p_dudt),
        copy(floes.p_dvdt),
        copy(floes.p_dξdt),
        copy(floes.p_dαdt),
        zeros(UInt8, nfloes),
    )
end

#=
Copy the fields that `timestep_floe_properties!` changes from `fixed_width_floes` back into
`floes`. `fixed_width_floes` must be on the CPU (use `adapt(Array, fixed_width_floes)`).
Fields are written through the StructArray's columns: iterating over `floes` would create
copies of each floe, so changes would be lost.
=#
function update_floes!(floes::StructArray{<:Floe{FT}}, fixed_width_floes::FixedWidthFloes{FT}) where FT
    fwf = fixed_width_floes
    for i in eachindex(floes)
        points = [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:fwf.n_points[i]]
        floes.poly[i] = GI.Polygon([GI.LinearRing(points)])::Polys{FT}
        floes.centroid[i] .= @view fwf.centroid[i, :]
        floes.stress_accum[i] .= @view fwf.stress_accum[i, :, :]
        floes.stress_instant[i] .= @view fwf.stress_instant[i, :, :]
        floes.strain[i] .= @view fwf.strain[i, :, :]
    end
    floes.height .= fwf.height
    floes.mass .= fwf.mass
    floes.moment .= fwf.moment
    floes.α .= fwf.α
    floes.u .= fwf.u
    floes.v .= fwf.v
    floes.ξ .= fwf.ξ
    floes.p_dxdt .= fwf.p_dxdt
    floes.p_dydt .= fwf.p_dydt
    floes.p_dudt .= fwf.p_dudt
    floes.p_dvdt .= fwf.p_dvdt
    floes.p_dξdt .= fwf.p_dξdt
    floes.p_dαdt .= fwf.p_dαdt
    return
end

"""
    Floe{FT}

Each model has a list of floes that are advected and updated throughout the simulation. The 
`Floe` struct holds the fields needed to compute all floe calculations. However, to create a floe,
only 2 arguments are neccessary and there are 2 more optional keyword arguments - all other fields
will be set to good starting values for a simualation.
 
!!! warning
    It is NOT reccomended for users to create individual floes using the `Floe` constructor. They should use the
    [`initialize_floe_field`](@ref) function to generate a field of floes and ensure that all of the fields are set
    correctly within the simulation. This function should only be used by developers who need to create an individual
    floe during the simulation (i.e. for fracturing and packing).

For those developers, here is how to construct a `Floe`:

    Floe{FT}(shape, height; floe_settings = FloeSettings(), rng = Xoshiro(), kwargs...)

The `FT` type should be provided as shown in the examples below:

## _Positional arguments_
- `shape::Union{<:Polys, <:PolyVec}`: either a polygon of type `Polys` or a vector of coordiantes of type `PolyVec` used to represent the shape of a floe.
- `height:FT`: height of the new floe to be created

## _Keyword arguments_
- `floe_settings::FloeSettings`: floe settings for a simulation (see [`FloeSettings`](@ref))
- `rng:<AbstractRNG`: random number generator to use to generate the subfloe points (if needed)
- `kwargs...`: additional field values can be passed in if desired

## _Examples_

- Creating a `Floe` with coordiantes:
```jldoctest floe
julia> coords = [[[0.0, 0.0], [0.0, 5.0], [10.0, 5.0], [10.0, 0.0], [0.0, 0.0]]];

julia> height = 1.0;

julia> Floe{Float64}(coords, height)
Floe{Float64}
  ⊢Centroid of (5.0, 2.5) m
  ⊢Height of 1.0 m
  ⊢Area of 50.0 m^2
  ∟Velocity of (u, v, ξ) of (0.0, 0.0, 0.0) in (m/s, m/s, rad/s)
```
- Creating a `Floe` with a polygon:
```jldoctest floe
julia> poly = make_polygon(coords, Float64);

julia> Floe{Float32}(poly, height; u = 1.0, ξ = 0.02)
Floe{Float32}
  ⊢Centroid of (5.0, 2.5) m
  ⊢Height of 1.0 m
  ⊢Area of 50.0 m^2
  ∟Velocity of (u, v, ξ) of (1.0, 0.0, 0.02) in (m/s, m/s, rad/s)
```

## _Fields_

!!! note
    There are quite a few fields in many different catagories, which will all be explained below.
    Eventually, it might be useful to make the `Floe` struct a "struct of structs" where each of these
    catagories is its own sub-struct. If you're interested in that, see the issue [#95](@ref).

The first catagory is **physical properties**. These have to do with the floe's physical shape. 

Before listing the fields, one important thing to know is that a floe's coordinates are represented by a `PolyVec`, which is a shorthand for a vector of a vector of a vector of floats. This sounds complicated, but it is simply a way of representing a polygon's coordinates. A Polygon's coordinates are of the form below, where the xy-coordinates are the exterior border of the floe and the wz-coordinates, or any other following sets of coordinates, describe holes within the floe:

```julia
coords = [
  [[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]],  # Exterior vertices of the polygon represented as a list of cartesian points
  [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]],  # Interior holes of the polygon represented as a list of cartesian points
  ...,  # Additional holes within the polygon represented as a list of cartesian points
 ]
 ```
 We will use the term `PolyVec` to describe this form of coordiantes and you will see it in the code if you take a look at the code base. It is also the form that floe coordinates are saved in output files.
 
| Physical Fields| Meaning                            | Type          |
| -------------- | ---------------------------------- | ------------- |
| poly           | polygon that represent's a floe's shape | Polys{Float64 or Float32}
| centroid       | floe's centroid                    | Float64 or Float32|
| height         | floe's height in [m]                 | Float64 or Float32|
| area           | floe's area in [m^2]                 | Float64 or Float32|
| mass           | floe's mass in [kg]                  | Float64 or Float32|
| rmax           | floe's maximum radius, the maximum <br> distance from centroid to vertex in [m] | Float64 or Float32|
| moment         | floe's mass moment of intertia in [kg m^2]    | Float64 or Float32|
| angles         | list of floe's vertex angles in [degrees] | Vector of Float64 or Float32|

The second catagory is **sub-floe points**. These are used for interpolation of the ocean and atmosphere onto the floe. They are a list of points within the floe. There are two ways they can be generated.
One is for them are randomly generated, with the user providing an initial target number of points. These are generated with a monte carlo generator (explained below).
The other way is for the points to be on a sub-grid within the floe. There are benefits and drawbacks to both strategies. 

| Sub-floe Point Fields| Meaning                                 | Type                        |
| ----------------- | --------------------------------------- | --------------------------- |
| x_subfloe_points   | floe's sub-floe points x-coordinates | Vector of Float64 or Float32|
| y_subfloe_points   | floe's sub-floe points points y-coordinates | Vector of Float64 or Float32|

The third catagory is **velocities and orientations**. Floe's have both linear and angular velocity and keep track of the angle that they have rotated since the begining of the simulation.
| Movement Fields| Meaning                         | Type               |
| -------------- | ------------------------------- | ------------------ |
| u              | floe's x-velocity in [m/s]        | Float64 or Float32 |
| v              | floe's x-velocity in [m/s]        | Float64 or Float32 |
| ξ              | floe's angular velocity in [rad/s]| Float64 or Float32 |
| α              | rotation from starting position<br> in [rad]| Float64 or Float32 |

The fourth catagory is **status**. These fields hold logistical information about each floe and through which process it originated.
| Status Fields  | Meaning                         | Type               |
| -------------- | ------------------------------- | ------------------ |
| status         | if the floe is still active in the simulation        | Subzero.Status (see below)|
| id             | unique floe id for tracking the floe throughout the simulation | Int |
| ghost_id       | if floe is not a ghost, `ghost_id = 0`, else it is in `[1, 4]`<br> as each floe can have up to 4 ghosts| Int |
| parent_id    | if floe is created from a fracture or the fusion of two floes, `parent_id` is a list of <br> the original floes' `id`, else it is emtpy | Vector of Ints |
| ghosts         | indices of floe's ghost floes within the floe list| Vector of Ints |

A Status object has two fields: a tag and a fuse index list. There are currently three different tags: `active`, `remove`, and `fuse`.
If a floe is `active`, it will continue in the simulation at the end of a timestep. If a floe's tag is `remove`, it will be removed at the end of the timestep.
This ususally happens if a floe exits the domain, or becomes unstable for some reason. If a floe is marked at `fuse`, this means that is is overlapping with another
floe by more than the user defined maximum overlap percent (see [`CollisionSettings`](@ref) for more information on this maximum overlap value.
If a floe is marked for fusion, the index of the floe it is supposed to fuse with will be listed in the `fuse_idx` list.  

The fifth catagory is **forces and collisions**. These fields hold information about the forces on each floe and the collisions it has been in.
| Force Fields      | Meaning                                    | Type                        |
| ----------------- | ------------------------------------------ | --------------------------- |
| fxOA              | x-force on floe from ocean and atmosphere in [N] | Float64 or Float32|
| fyOA              | y-force on floe from ocean and atmosphere in [N] | Float64 or Float32|
| trqOA             | torque on floe from ocean and atmosphere in [N m]| Float64 or Float32|
| hflx_factor       | coefficent of floe height to get heat flux directly <br> under floe in [W/m^3]| Float64 or Float32|
| overarea          | total overlap of floe from collisions in [m^2]   | Float64 or Float32|
| collision_force   | forces on floe from collisions in [N]            | Float64 or Float32|
| collision_trq     | torque on floe from collisions in [N m]          | Float64 or Float32|
| interactions      | each row holds one collision's information, see below for more information | `n`x7 Matrix of Float64 or Float32 <br> where `n` is the number of collisions|
| stress_accum      | stress accumulated over the floe over past timesteps given StressCalculator, where it is of the form [xx yx; xy yy] | 2x2 Matrix{AbstractFloat}|
| stress_instant    | instantaneous stress on floe in current timestep | 2x2 Matrix{AbstractFloat} 
| strain            | strain on floe where it is of the form [ux vx; uy vy] | 2x2 Matrix of Float64 or Float32|
| damage             | damage a floe assumulates through interactions (currently unused; see [`DamageStressCalculator`](@ref))|  Float64 or Float32|

The `interactions` field is a matrix where every row is a different collision that the floe has experienced in the given timestep. There are then seven columns, which are as follows:
- `floeidx`, which is the index of the floe that the current floe collided with
- `xforce`, which is the force in the x-direction caused by the collision
- `yforce`, which is the force in the y-direction caused by the collision
- `xpoint`, which is the x-coordinate of the collision point, which is the x-centroid of the overlap between the floes
- `ypoint`, which is the y-coordinate of the collision point, which is the y-centroid of the overlap between the floes
- `torque`, which is the torque caused by the collision
- `overlap`, which is the overlap area between the two floes in the collision. 
You can use these column names to access columns of `interactions`. For example: `floe.interactions[:, xforce]` gives the list of x-force values for all collisions a given floe was involved in as that timestep.

The fifth catagory is **timesteping values**. These are values from the previous timestep that are needed for the current timestepping scheme.
| Previous Value Fields | Meaning                                       | Type               |
| --------------------- | --------------------------------------------- | -------------------|
| p_dxdt                | previous timestep x-velocity (u) in [m/s]         | Float64 or Float32|
| p_dydt                | previous timestep y-velocity (v) in [m/s]         | Float64 or Float32|
| p_dudt                | previous timestep x-acceleration in [m/s^2]   | Float64 or Float32|
| p_dvdt                | previous timestep y-acceleration in [m/s^2]   | Float64 or Float32|
| p_dαdt                | previous timestep angular-velocity in [rad/s] | Float64 or Float32|
| p_dξdt                | previous timestep time angular acceleration in [rad/s^2] | Float64 or Float32|
"""
function Floe{FT}(shape, height; floe_settings = FloeSettings(), rng = Xoshiro(), kwargs...) where FT
    poly = _correct_floe_poly(FT, shape)
    # Floe physical properties
    centroid = collect(centroid_poly(poly, FT))
    area = area_poly(poly, FT)
    mass = area * height * floe_settings.ρi
    moment = _calc_moment_inertia(FT, poly, centroid, height; ρi = floe_settings.ρi)
    rmax = _calc_max_radius(poly, centroid, FT)
    angles = angles_poly(poly, FT)
    # Generate Monte Carlo points
    status = Status()
    x_subfloe_points, y_subfloe_points, status = generate_subfloe_points(floe_settings.subfloe_point_generator, poly, centroid, area, status, rng)
    # Generate status
    return Floe{FT}(; poly, height, centroid, area, mass, rmax, moment, angles, status, x_subfloe_points, y_subfloe_points, kwargs...)
end

# if no type is provided, it will default to Float64 - not a argument as this should not be user facing!
Floe(args...; kwargs...) = Floe{Float64}(args...; kwargs...)

# Pretty printing for Floe showing key physical fields
function Base.show(io::IO, floe::Floe{FT}; digits = 5) where FT
    overall_summary = "Floe{$FT}"

    floe_centroid_summary = "Centroid of ($(round(floe.centroid[1]; digits)), $(round(floe.centroid[2]; digits))) m"
    floe_height_summary = "Height of $(round(floe.height; digits)) m"
    floe_area_summary = "Area of $(round(floe.area; digits)) m^2"
    floe_velocity_summary = "Velocity of (u, v, ξ) of ($(round(floe.u; digits)), $(round(floe.v; digits)), $(round(floe.ξ; digits))) in (m/s, m/s, rad/s)"
    
    print(io, overall_summary, "\n",
        "  ⊢", floe_centroid_summary, "\n",
        "  ⊢", floe_height_summary, "\n",
        "  ⊢", floe_area_summary, "\n",
        "  ∟", floe_velocity_summary, "\n")
end
