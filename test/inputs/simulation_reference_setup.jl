#=
The simulation of test_simulation_reference.jl, shared by that test and by
make_simulation_reference.jl, which makes its reference results with the original CPU
implementation (Subzero at ae96ff2, before the GPU port). So everything here must also work
with that version of Subzero.

It's a small version of the shear flow example: periodic boundaries (so ghost floes),
collisions, and coupling with a shear ocean current. The floes are saved as plain arrays,
because a Voronoi floe field isn't reproducible, and because any Subzero version can load
plain arrays.
=#

# Floe fields that are saved and loaded as they are
const REFERENCE_FLOE_FIELDS = (
    :centroid, :height, :area, :mass, :rmax, :moment, :angles,
    :α, :u, :v, :ξ, :fxOA, :fyOA, :trqOA, :hflx_factor, :overarea,
    :collision_force, :collision_trq, :interactions, :num_inters,
    :stress_accum, :stress_instant, :strain, :damage,
    :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt,
    :x_subfloe_points, :y_subfloe_points, :id,
)

# Floe fields that are compared with the reference results
const REFERENCE_RESULT_FIELDS = (
    :id, :centroid, :height, :mass, :moment, :α, :u, :v, :ξ,
    :stress_accum, :stress_instant, :strain,
)

# Timesteps after which the floes are compared with the reference results
const REFERENCE_TIMESTEPS = (1, 10, 50, 100)

const REFERENCE_GRID = (; x0 = 0.0, xf = 2e4, y0 = 0.0, yf = 2e4, Δx = 2e3, Δy = 2e3)

reference_grid() = RegRectilinearGrid(; REFERENCE_GRID...)

reference_domain(grid) = Domain(;
    north = PeriodicBoundary(North; grid), south = PeriodicBoundary(South; grid),
    east = PeriodicBoundary(East; grid), west = PeriodicBoundary(West; grid),
)

reference_floe_settings(grid) =
    FloeSettings(subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2))

# The floes as a Dict of plain arrays: the exterior ring of each floe as a vector of [x, y],
# and the fields in REFERENCE_FLOE_FIELDS
function reference_floes_to_dict(floes)
    d = Dict{String, Any}(string(f) => collect(getproperty(floes, f)) for f in REFERENCE_FLOE_FIELDS)
    d["coords"] = [[collect(GI.getcoord(p)) for p in GI.getpoint(GI.getexterior(poly))]
        for poly in floes.poly]
    return d
end

# The floes saved by reference_floes_to_dict
reference_floes_from_dict(d) = StructArray([
    Floe{Float64}(;
        poly = Subzero.make_polygon([[copy(c) for c in d["coords"][i]]]),
        (f => deepcopy(d[string(f)][i]) for f in REFERENCE_FLOE_FIELDS)...,
    ) for i in eachindex(d["coords"])
])

#=
The simulation with the given floes. `kwargs` are passed on to `Simulation`, e.g. the
backend, which the original implementation doesn't have.
=#
function reference_simulation(floes; kwargs...)
    grid = reference_grid()
    domain = reference_domain(grid)
    half_nx = fld(grid.Nx, 2)
    u_vec = [range(0, 0.5, length = half_nx); 0.5; range(0.5, 0, length = half_nx)]
    uvels = repeat(transpose(u_vec), outer = (grid.Ny + 1, 1))
    ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)
    atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)
    model = Model(; grid, ocean, atmos, domain, floes)
    consts = Constants(; E = 1.5e3 * (mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area))))
    return Simulation(; model, consts, Δt = 20, nΔt = maximum(REFERENCE_TIMESTEPS),
        floe_settings = reference_floe_settings(grid), rng = Xoshiro(1), kwargs...)
end

#=
Run the simulation and return the floes after each timestep in REFERENCE_TIMESTEPS, as a
Dict from the timestep to a Dict of the fields in REFERENCE_RESULT_FIELDS and the floes'
coordinates. Logging is switched off.
=#
function run_reference_simulation(sim)
    results = Dict{Int, Dict{String, Any}}()
    with_logger(NullLogger()) do
        for tstep in 0:sim.nΔt
            timestep_sim!(sim, tstep)
            if tstep in REFERENCE_TIMESTEPS
                floes = sim.model.floes
                r = Dict{String, Any}(string(f) => deepcopy(collect(getproperty(floes, f)))
                    for f in REFERENCE_RESULT_FIELDS)
                r["coords"] = [[collect(GI.getcoord(p)) for p in GI.getpoint(GI.getexterior(poly))]
                    for poly in floes.poly]
                results[tstep] = r
            end
        end
    end
    return results
end
