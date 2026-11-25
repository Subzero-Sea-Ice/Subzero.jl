using Subzero, StructArrays, Statistics, JLD2
import LibGEOS as LG

# User Inputs
const FT = Float64

const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 20
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = RegRectilinearGrid(FT; x0 = -Lx, xf = Lx, y0 = -Ly, yf = Ly, Δx = Δgrid, Δy = Δgrid)

ocean = Ocean(FT; grid, u = -0.2, v = 0.0, temp = -1.0)

atmos = Atmos(FT; grid, u = 0.0, v = 0.0, temp = -3.0)

# Domain creation - boundaries and topography
nboundary = OpenBoundary(North, FT; grid)
sboundary = OpenBoundary(South, FT; grid)
eboundary = OpenBoundary(East, FT; grid)
wboundary = OpenBoundary(West, FT; grid)

domain = Subzero.Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# Floe instantiation
nfloes = 100
file = jldopen("examples/floe_shapes.jld2", "r")
nfloes = nfloes > size(file["floe_vertices"], 1) ? size(file["floe_vertices"], 1) : nfloes
floe_coords = file["floe_vertices"][1:1]
floe_arr = initialize_floe_field(FT, floe_coords, domain, hmean, Δh)
close(file)

model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(FT; E = modulus)
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(; model, consts, Δt, nΔt = 4000, verbose = true)

# Output setup
# initwriter = InitialStateOutputWriter(dir = "output/voronoi", filename = "initial_state.jld2", overwrite = true)
# floewriter = FloeOutputWriter(10, dir = "output/voronoi", filename = "f.jld2", overwrite = true)

# # Run simulation
# run!(simulation, [initwriter, floewriter])