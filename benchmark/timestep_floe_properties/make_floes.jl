# Generate the benchmark input: the docs example `shear_flow.jl` on a larger domain with
# more floes, spun up for NSPINUP timesteps so velocities, forces, and interactions are
# non-trivial. Floe fields are saved as plain arrays so that any Subzero version can load
# them (a Voronoi floe field isn't reproducible, so generate it once and reuse the file).
#
# Usage: julia --project=env-branch make_floes.jl [nfloes] [nspinup] [outfile]
using Subzero, JLD2, Random, Statistics
import GeoInterface as GI

const NFLOES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
const NSPINUP = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100
const OUTFILE = length(ARGS) >= 3 ? ARGS[3] : "floes.jld2"

const FT = Float64
# Same floe density as shear_flow.jl (50 floes in 1e5 × 1e5 m), so floes have the same size
# (rounded to an even number of grid cells, which the shear profile below needs)
const Δgrid = 2e3
const Lx = 2Δgrid * round(1e5 * sqrt(NFLOES / 50) / 2Δgrid)
const Ly = Lx
const hmean = 0.25
const Δh = 0.0
const Δt = 20

grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
domain = Domain(;
    north = PeriodicBoundary(North; grid), south = PeriodicBoundary(South; grid),
    east = PeriodicBoundary(East; grid), west = PeriodicBoundary(West; grid),
)
half_nx = fld(grid.Nx, 2)
u_vec = [range(0, 0.5, length = half_nx); 0.5; range(0.5, 0, length = half_nx)]
uvels = repeat(transpose(u_vec), outer = (grid.Ny + 1, 1))
ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)
atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)

floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2))
floe_generator = VoronoiTesselationFieldGenerator(; nfloes = NFLOES, concentrations = [0.75], hmean, Δh)
@info "Generating floe field" NFLOES Lx
@time floe_arr = initialize_floe_field(FT; generator = floe_generator, domain, rng = Xoshiro(1), floe_settings)
@info "Generated floes" length(floe_arr)

model = Model(; grid, ocean, atmos, domain, floes = floe_arr)
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(; E = modulus)
simulation = Simulation(; model, consts, Δt, nΔt = NSPINUP, floe_settings, rng = Xoshiro(1))
@info "Spinning up" NSPINUP
@time run!(simulation)

floes = simulation.model.floes
field(f) = collect(getproperty(floes, f))
jldsave(OUTFILE;
    Δt, Lx, Ly, nspinup = NSPINUP,
    # exterior ring of each floe as a vector of [x, y]
    coords = [[collect(GI.getcoord(p)) for p in GI.getpoint(GI.getexterior(poly))] for poly in floes.poly],
    (f => field(f) for f in (
        :centroid, :height, :area, :mass, :rmax, :moment, :angles,
        :α, :u, :v, :ξ, :fxOA, :fyOA, :trqOA, :hflx_factor, :overarea,
        :collision_force, :collision_trq, :interactions, :num_inters,
        :stress_accum, :stress_instant, :strain, :damage,
        :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt,
    ))...,
)
@info "Saved" OUTFILE length(floes) mean(floes.num_inters) maximum(floes.num_inters) mean(GI.npoint.(floes.poly))
