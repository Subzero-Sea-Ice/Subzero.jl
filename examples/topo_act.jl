using JLD2, Random, Statistics, Subzero, CairoMakie, GeoInterfaceMakie
using NCDatasets, Interpolations, StaticArrays, LinearAlgebra
using JSON
import GeometryBasics as GB

# User Inputs
const FT = Float64
const Lx = 512e3
const Ly = 512e3
const Δgrid = 2e3
const hmean = 1.78
const Δh = 0.0
const Δt = 1.0

# Model initialization
grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
x = range(grid.x0, stop=grid.xf, length=grid.Nx + 1)
y = range(grid.y0, stop=grid.yf, length=grid.Ny + 1)


# Reading wind data
ncfile_wind = "/home/pvishwaja/poly_topo/APRIL_2016_Wind_box_2km_all_times.nc"
ds = NCDataset(ncfile_wind, "r")
xatm = ds["x"][:] .* 1e3
yatm = ds["y"][:] .* 1e3
uatm = ds["u"][:, :, :]   #(x, y, time)
vatm = ds["v"][:, :, :]
uatm_perm = permutedims(uatm, (3, 1, 2))  # (time, x, y)
vatm_perm = permutedims(vatm, (3, 1, 2))
time_atm = ds["time"][:] ./ 8.0
close(ds)
itp_uatm = interpolate((time_atm, xatm, yatm),uatm_perm,Gridded(Linear()))
itp_vatm = interpolate((time_atm, xatm, yatm),vatm_perm,Gridded(Linear()))
t0_days = first(time_atm)
uatm_init = Array{Float64}(undef, length(x), length(y))
vatm_init = Array{Float64}(undef, length(x), length(y))
for i in 1:length(x)
    for j in 1:length(y)
        uatm_init[i, j] = itp_uatm(t0_days, x[i], y[j])
        vatm_init[i, j] = itp_vatm(t0_days, x[i], y[j])
    end
end
println("Read ocean data: size(u) = ", size(uatm), "First time = ", t0_days)


# Reading ocean data
ncfile_oce = "/home/pvishwaja/poly_topo/APRIL_2016_Ocean_box_2km_all_times.nc"
ds = NCDataset(ncfile_oce, "r")
xoc = ds["x"][:] .* 1e3
yoc = ds["y"][:] .* 1e3
time_oce = Float64.(ds["time"][:])     
uoce = ds["u"][:, :, :]     
voce = ds["v"][:, :, :]
uoce_perm = permutedims(uoce, (3, 1, 2))  # (time, x, y)
voce_perm = permutedims(voce, (3, 1, 2))
close(ds)
itp_uoce = interpolate((time_oce, xoc, yoc),uoce_perm,Gridded(Linear()))
itp_voce = interpolate((time_oce, xoc, yoc),voce_perm,Gridded(Linear()))
t0_days = first(time_oce)
uoce_init = Array{Float64}(undef, length(x), length(y))
voce_init = Array{Float64}(undef, length(x), length(y))
for i in 1:length(x)
    for j in 1:length(y)
        uoce_init[i, j] = itp_uoce(t0_days, x[i], y[j])
        voce_init[i, j] = itp_voce(t0_days, x[i], y[j])
    end
end
println("Read ocean data: size(u) = ", size(uoce), "First time = ", t0_days)


#Read uvels, vvels from netCDF file
atmos = Atmos(; grid, u = uatm_init, v = vatm_init, temp = -1.0)
ocean = Ocean(; grid, u = uoce_init, v = uoce_init, temp = 0.0)

# Domain creation
nboundary = OpenBoundary(North; grid)
sboundary  = OpenBoundary(South; grid)
eboundary = OpenBoundary(East; grid)
#wboundary = OpenBoundary(West; grid)
wboundary = MovingBoundary(West; grid, u = 0.0, v = 0.0)

#topo1_raw = JSON.parsefile("/home/pvishwaja/poly_topo/topo2_25m.json")
#topo1 = Vector{Vector{Vector{Float64}}}(topo1_raw)
#topo2_raw = JSON.parsefile("/home/pvishwaja/poly_topo/topo2_25m.json")
#topo2 = Vector{Vector{Vector{Float64}}}(topo2_raw)
#topo1_km = [ [ p .* 1000 for p in poly ] for poly in topo1 ]
#topo2_km = [ [ p .* 1000 for p in poly ] for poly in topo2 ]
#coords1 = [topo1_km,topo2_km]

topo18_raw = JSON.parsefile("/home/pvishwaja/poly_topo/topo18.json")
topo18 = Vector{Vector{Vector{Float64}}}(topo18_raw)
topo18_km = [[p .* 1000 for p in topo] for topo in topo18]
coords = [[topo] for topo in topo18_km]


topo_arr = initialize_topography_field(FT; coords)

#domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary, topography = topo_arr)
coupling_settings = CouplingSettings(two_way_coupling_on = false,)
bond_settings = BondSettings(bonds_on = true, k_b = 1e6, c_b = 1e8, λ_Rb = 1.0, λ_lb = 0.8, SF_b = 0.077, ζ_b = 0.7,)
bond_fracture_settings = Bond_FractureSettings(σt_max_bond = 30000, σc_max_bond = 1e6, τ_max_bond = 1e6,)
floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(grid, 2),)

floe_arr = initialize_floe_field(
    FT,
    2500,
    [1.0],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

function polygon_edges(poly_coords::Vector{Vector{Float64}})
    n = length(poly_coords)
    pts = [GB.Point2(c[1], c[2]) for c in poly_coords]
    [GB.Line(pts[i], pts[mod1(i+1, n)]) for i in 1:n]
end

function edges_equal(e1::GB.Line, e2::GB.Line; rtol=1e-3)
    scale = 0.5 * (norm(e1.points[2] - e1.points[1]) + norm(e2.points[2] - e2.points[1]))
    atol = rtol * scale  # dynamic absolute tolerance

    (isapprox(e1.points[1], e2.points[1]; atol=atol, rtol=rtol) &&
     isapprox(e1.points[2], e2.points[2]; atol=atol, rtol=rtol)) ||
    (isapprox(e1.points[1], e2.points[2]; atol=atol, rtol=rtol) &&
     isapprox(e1.points[2], e2.points[1]; atol=atol, rtol=rtol))
end

struct BoundingBox
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end
num_floes = length(floe_arr.coords)
function compute_bbox(poly_coords::Vector{Vector{Float64}})
    xs = [c[1] for c in poly_coords]
    ys = [c[2] for c in poly_coords]
    BoundingBox(minimum(xs), maximum(xs), minimum(ys), maximum(ys))
end

bboxes = [compute_bbox(floe_arr.coords[i][1]) for i in 1:num_floes]
function bbox_intersect(b1::BoundingBox, b2::BoundingBox; atol=1e-8)
    return !(
        (b1.xmax <= b2.xmin - atol) ||  (b1.xmin >= b2.xmax + atol) ||
        (b1.ymax <= b2.ymin - atol) ||  (b1.ymin >= b2.ymax + atol))
end

bonds = Vector{Vector{Int}}(undef, num_floes)

for i_idx in eachindex(floe_arr)
    edges_i = polygon_edges(floe_arr.coords[i_idx][1])
    num_edges = length(edges_i)
    bonds[i_idx] = zeros(Int, num_edges - 1)
end

for i_idx in eachindex(floe_arr)
    floe = floe_arr[i_idx]
    coords_rearr = [c for c in floe_arr.coords]
    edges_i = polygon_edges(coords_rearr[i_idx][1])
    num_edges_i = length(edges_i) - 1
    for j_idx in i_idx+1:num_floes
        if !bbox_intersect(bboxes[i_idx], bboxes[j_idx])
            continue
        end
        edges_j = polygon_edges(coords_rearr[j_idx][1])
        num_edges_j = length(edges_j) - 1

        for (e_idx_i, e_i) in enumerate(edges_i[1:num_edges_i])
            if bonds[i_idx][e_idx_i] != 0
                continue
            end
            for (e_idx_j, e_j) in enumerate(edges_j[1:num_edges_j])
                if bonds[j_idx][e_idx_j] != 0
                    continue
                end

                if edges_equal(e_i, e_j)
                    bonds[i_idx][e_idx_i] = floe_arr.id[j_idx]
                    bonds[j_idx][e_idx_j] = floe_arr.id[i_idx]
                    break
                end
            end
        end
    end
    floe.bonds= bonds[i_idx]
    floe_arr[i_idx] = floe
end

function point_in_poly(x, y, xs, ys)
    inside = false
    j = length(xs)
    for i in eachindex(xs)
        if ((ys[i] > y) != (ys[j] > y)) &&
           (x < (xs[j] - xs[i]) * (y - ys[i]) / (ys[j] - ys[i] + eps()) + xs[i])
            inside = !inside
        end
        j = i
    end
    return inside
end

ncfile_thick = "/home/pvishwaja/poly_topo/APRIL_2016_Thick_box_2km.nc"
ds = NCDataset(ncfile_thick)
sithick = ds["sithick"][:, :, 1]   
Xt = repeat(x', length(y), 1)   
Yt = repeat(y, 1, length(x))   
for i_idx in eachindex(floe_arr)
    floe = floe_arr[i_idx]
    verts = floe.coords[1]
    xs = [v[1] for v in verts]   # x-coordinates
    ys = [v[2] for v in verts]   # y-coordinates
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    mask_bbox = (Xt .>= xmin) .& (Xt .<= xmax) .& (Yt .>= ymin) .& (Yt .<= ymax)
    idxs = findall(mask_bbox)
    thick_vals = Float64[]
    for idx in idxs
        xi = Xt[idx]
        yi = Yt[idx]
        if point_in_poly(xi, yi, xs, ys)
            val = sithick[idx]
            if !isnan(val)
                push!(thick_vals, val)
            end
        end
    end
    floe.height = isempty(thick_vals) ? 0.0 : mean(thick_vals)
    floe_arr[i_idx] = floe
end
function smooth_floe_thickness!(floe_arr; β=0.2, niter=3)
    hnew = similar(floe_arr, Float64)

    for iter in 1:niter
        for i in eachindex(floe_arr)
            hi = floe_arr.height[i]
	    nbrs = [j for j in floe_arr.bonds[i] if j > 0]
            if isempty(nbrs)
                hnew[i] = hi
            else
                h_mean = mean(floe_arr.height[j] for j in nbrs)
                hnew[i] = hi + β *(h_mean - hi)
            end
        end

        for i in eachindex(floe_arr)
            floe = floe_arr[i]
            floe.height = hnew[i]
            floe_arr[i] = floe
        end
    end
end
smooth_floe_thickness!(floe_arr)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
f = 1.46e-4
consts = Constants(E = modulus, f = f)

# Run simulation
run_time!(simulation) =  @time run!(simulation)
dir = "/scratch/pvishwaja/poly_topo/check_dt1"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(1, dir = dir, overwrite = true)
gridwriter = GridOutputWriter(1, grid, (grid.Nx, grid.Ny), dir = dir, overwrite = true)
checkpointwriter = CheckpointOutputWriter(1, dir = dir, filename = "checkpoint.jld2", overwrite = true)
writers = OutputWriters(initwriter, floewriter, checkpointwriter, gridwriter )
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 10,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
    bond_settings = bond_settings,
    bond_fracture_settings = bond_fracture_settings,
)

function run_with_change_atm_oce(simulation)
    Subzero.startup_sim(simulation, nothing, 1)
    uoce = similar(simulation.model.ocean.u)
    voce = similar(simulation.model.ocean.v)
    uatm = similar(simulation.model.atmos.u)
    vatm = similar(simulation.model.atmos.v)
    tstep = 0
    while tstep <= simulation.nΔt
        time_sec  = tstep * simulation.Δt
        time_days = time_sec / 86400.0
        for i in eachindex(x)
            for j in eachindex(y)
                uoce[i, j] = itp_uoce(time_days, x[i], y[j])
                voce[i, j] = itp_voce(time_days, x[i], y[j])
            end
        end
        simulation.model.ocean.u .= uoce
        simulation.model.ocean.v .= voce
        for i in eachindex(x)
            for j in eachindex(y)
                uatm[i, j] = itp_uatm(time_days, x[i], y[j])
                vatm[i, j] = itp_vatm(time_days, x[i], y[j])
            end
        end
        simulation.model.atmos.u .= uatm
        simulation.model.atmos.v .= vatm
        timestep_sim!(simulation, tstep)
        tstep += 1
    end
    Subzero.teardown_sim(simulation)
end

run_with_change_atm_oce(simulation)

plot_sim(
    "/scratch/pvishwaja/poly_topo/check_dt1/floes.jld2",
    "/scratch/pvishwaja/poly_topo/check_dt1/initial_state.jld2",
    Δt,
    "/scratch/pvishwaja/poly_topo/check_dt1/wind_flow.mp4")
