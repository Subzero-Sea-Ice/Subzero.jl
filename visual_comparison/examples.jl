#=
The three documentation examples that are compared: their setup, as in
docs/literate/examples/, and how their floes are saved and loaded.

This file is used with two versions of Subzero: the original CPU implementation and the
current one. Everything here must work with both, so it only uses the API the examples use.

The examples make their floes with a Voronoi tessellation, which isn't reproducible with a
fixed rng. So the floes are made once (make_floes.jl), saved as plain arrays, and both
versions start from the saved floes.
=#

const EXAMPLES = ("shear_flow", "simple_strait", "forcing_contained_floes", "moving_bounds")

# The titles of the examples in the documentation
const EXAMPLE_TITLES = Dict(
    "shear_flow" => "Shear Flow with Periodic Boundaries",
    "simple_strait" => "Simple Strait Simulation",
    "forcing_contained_floes" => "Floes Bounded by Ocean Currents",
    "moving_bounds" => "Moving Boundaries",
)

const FT = Float64
const Δt = 20

# Floe fields that are saved and loaded as they are, besides the polygon
const FLOE_FIELDS = (
    :centroid, :height, :area, :mass, :rmax, :moment, :angles,
    :α, :u, :v, :ξ, :fxOA, :fyOA, :trqOA, :hflx_factor, :overarea,
    :collision_force, :collision_trq, :interactions, :num_inters,
    :stress_accum, :stress_instant, :strain, :damage,
    :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt,
    :x_subfloe_points, :y_subfloe_points, :id,
)

#=
The setup of an example: the model without floes, how its floes are made (`init_floes!`
changes them after they're made), the keyword arguments of its `Constants` besides `E`, and
the settings of its simulation (a function of the floes, because the fracture criteria depend
on them).
=#
example_setup(name::AbstractString) = example_setup(Val(Symbol(name)))

grid_100km() = RegRectilinearGrid(; x0 = 0.0, xf = 1e5, y0 = 0.0, yf = 1e5, Δx = 2e3, Δy = 2e3)

# docs/literate/examples/shear_flow.jl
function example_setup(::Val{:shear_flow})
    grid = grid_100km()
    domain = Domain(;
        north = PeriodicBoundary(North; grid), south = PeriodicBoundary(South; grid),
        east = PeriodicBoundary(East; grid), west = PeriodicBoundary(West; grid),
    )
    half_nx = fld(grid.Nx, 2)
    u_vec = [range(0, 0.5, length = half_nx); 0.5; range(0.5, 0, length = half_nx)]
    uvels = repeat(transpose(u_vec), outer = (grid.Ny + 1, 1))
    ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)
    atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)
    floe_settings = FloeSettings(
        subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2))
    generator = VoronoiTesselationFieldGenerator(; nfloes = 50, concentrations = [0.75],
        hmean = 0.25, Δh = 0.0)
    return (; grid, domain, ocean, atmos, floe_settings, generator, rng_seed = 1,
        floe_kwargs = (;), init_floes! = floes -> nothing, consts_kwargs = (;),
        settings = floes -> (;))
end

# docs/literate/examples/simple_strait.jl
function example_setup(::Val{:simple_strait})
    grid = grid_100km()
    island1 = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
    island2 = [[[2e4, 6e4], [4e4, 6.5e4], [4.5e4, 6.5e4], [4.5e4, 6e4], [2e4, 6e4]]]
    topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
    topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]
    topography = initialize_topography_field(FT; coords = [island1, island2, topo1, topo2])
    domain = Domain(;
        north = PeriodicBoundary(North; grid), south = PeriodicBoundary(South; grid),
        east = CollisionBoundary(East; grid), west = CollisionBoundary(West; grid),
        topography,
    )
    ocean = Ocean(; grid, u = 0.0, v = -0.3, temp = 0.0)
    atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)
    floe_settings = FloeSettings(
        subfloe_point_generator = MonteCarloPointsGenerator(; npoints = 100, ntries = 10, err = 0.1),
        stress_calculator = DecayAreaScaledCalculator(),
    )
    generator = VoronoiTesselationFieldGenerator(; nfloes = 75, concentrations = [0.7],
        hmean = 0.25, Δh = 0.0)
    settings(floes) = (;
        fracture_settings = FractureSettings(; fractures_on = true,
            criteria = HiblerYieldCurve(floes), Δt = 75, npieces = 3, deform_on = false),
        ridgeraft_settings = RidgeRaftSettings(; ridge_raft_on = true, Δt = 10),
    )
    return (; grid, domain, ocean, atmos, floe_settings, generator, rng_seed = 3,
        floe_kwargs = (;), init_floes! = floes -> nothing, consts_kwargs = (;), settings)
end

# docs/literate/examples/forcing_contained_floes.jl
function example_setup(::Val{:forcing_contained_floes})
    grid = grid_100km()
    domain = Domain(;
        north = OpenBoundary(North; grid), south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid), west = OpenBoundary(West; grid),
    )
    ocean_uvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))
    ocean_vvels = zeros(FT, (grid.Nx + 1, grid.Ny + 1))
    for i in CartesianIndices(ocean_vvels)
        r, c = Tuple(i)
        if r ≤ 5
            ocean_uvels[i] = 0.6
        elseif r ≥ grid.Nx - 4
            ocean_uvels[i] = -0.6
        end
        if c ≤ 5
            ocean_vvels[i] = 0.6
        elseif c ≥ grid.Ny - 4
            ocean_vvels[i] = -0.6
        end
    end
    ocean_uvels[10:40, 20:30] .= 0.3
    ocean = Ocean(; grid, u = ocean_uvels, v = ocean_vvels, temp = 0)
    atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)
    floe_settings = FloeSettings(;
        subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2))
    generator = VoronoiTesselationFieldGenerator(; nfloes = 300, concentrations = [0.4],
        hmean = 0.25, Δh = 0.0)
    Lx = Ly = 1e5
    floe_bounds = Subzero.make_polygon([[[0.1Lx, 0.1Ly], [0.1Lx, 0.9Ly], [0.9Lx, 0.9Ly],
        [0.9Lx, 0.1Ly], [0.1Lx, 0.1Ly]]], FT)
    return (; grid, domain, ocean, atmos, floe_settings, generator, rng_seed = 1,
        floe_kwargs = (; floe_bounds), init_floes! = floes -> nothing, consts_kwargs = (;),
        settings = floes -> (;))
end

#= docs/literate/examples/moving_bounds.jl. The example defines weld settings but doesn't pass
them to its simulation, so there's no welding here either. =#
function example_setup(::Val{:moving_bounds})
    grid = grid_100km()
    domain = Domain(;
        north = MovingBoundary(North; grid, u = 0.0, v = -0.1),
        south = MovingBoundary(South; grid, u = 0.0, v = 0.1),
        east = PeriodicBoundary(East; grid), west = PeriodicBoundary(West; grid),
    )
    ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
    atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = -1.0)
    floe_settings = FloeSettings()
    generator = VoronoiTesselationFieldGenerator(; nfloes = 100, concentrations = [1.0],
        hmean = 0.25, Δh = 0.125)
    function init_floes!(floes)
        floes.u .= 0
        floes.v .= -0.01
        return
    end
    settings(floes) = (; ridgeraft_settings = RidgeRaftSettings(; ridge_raft_on = true,
        Δt = 150, domain_gain_probability = 0.5))
    return (; grid, domain, ocean, atmos, floe_settings, generator, rng_seed = 1,
        floe_kwargs = (;), init_floes!, consts_kwargs = (; Cd_io = 0.0, f = 0.0, turnθ = 0.0),
        settings)
end

# The example's floes, made as in the example
function make_example_floes(name)
    s = example_setup(name)
    floes = initialize_floe_field(FT; generator = s.generator, domain = s.domain,
        rng = Xoshiro(s.rng_seed), floe_settings = s.floe_settings, s.floe_kwargs...)
    s.init_floes!(floes)
    return floes
end

# The floes as a Dict of plain arrays: the exterior ring of each floe as a vector of [x, y],
# and the fields in FLOE_FIELDS
function floes_to_dict(floes)
    d = Dict{String, Any}(string(f) => collect(getproperty(floes, f)) for f in FLOE_FIELDS)
    d["coords"] = [[collect(GI.getcoord(p)) for p in GI.getpoint(GI.getexterior(poly))]
        for poly in floes.poly]
    return d
end

# The floes saved by floes_to_dict
floes_from_dict(d) = StructArray([
    Floe{FT}(;
        poly = Subzero.make_polygon([[copy(c) for c in d["coords"][i]]]),
        (f => deepcopy(d[string(f)][i]) for f in FLOE_FIELDS)...,
    ) for i in eachindex(d["coords"])
])

#=
The example's simulation, with the given floes, running for `nsteps` timesteps. `kwargs`
are passed on to `Simulation`, e.g. the backend, which the original implementation
doesn't have.
=#
function example_simulation(name, floes, nsteps; kwargs...)
    s = example_setup(name)
    model = Model(; grid = s.grid, ocean = s.ocean, atmos = s.atmos, domain = s.domain, floes)
    consts = Constants(; E = 1.5e3 * (mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area))),
        s.consts_kwargs...)
    return Simulation(; model, consts, Δt, nΔt = nsteps, floe_settings = s.floe_settings,
        rng = Xoshiro(1), s.settings(floes)..., kwargs...)
end

# The exterior rings of polygons as vectors of (x, y)
ring_coords(polys) = [[(Float64(GI.x(p)), Float64(GI.y(p))) for p in GI.getpoint(GI.getexterior(poly))]
    for poly in polys]
