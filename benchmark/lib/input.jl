# Load the floes saved by make_input.jl, and make the model they live in. This works with
# any Subzero version, so the baseline and the current version get exactly the same input.

const INFILE = get(ENV, "BENCH_FLOES", "floes_10k.jld2")

# Fields saved by make_input.jl that are copied into the floes as they are
const SAVED_FIELDS = (
    :height, :area, :mass, :rmax, :moment, :angles,
    :α, :u, :v, :ξ, :fxOA, :fyOA, :trqOA, :hflx_factor, :overarea,
    :collision_force, :collision_trq, :num_inters,
    :stress_accum, :stress_instant, :strain, :damage,
    :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt,
)

const INPUT = Ref{Dict{String, Any}}()
function input()
    isassigned(INPUT) && return INPUT[]
    d = load(INFILE)
    # Input files made before make_input.jl saved them have no grid spacing or subfloe
    # points. The subfloe points of a SubGridPointsGenerator don't depend on the rng, so they
    # can be made again from the polygons.
    get!(d, "Δgrid", 2e3)
    get!(d, "npoint_per_cell", 2)
    if !haskey(d, "x_subfloe_points")
        generator = SubGridPointsGenerator(Float64; Δg = d["Δgrid"] / d["npoint_per_cell"])
        points = map(d["coords"], d["centroid"], d["area"]) do coords, centroid, area
            poly = make_polygon([[collect(c) for c in coords]])
            x, y, _ = Subzero.generate_subfloe_points(generator, poly, centroid, area, Subzero.Status(), Xoshiro(1))
            (x, y)
        end
        d["x_subfloe_points"] = first.(points)
        d["y_subfloe_points"] = last.(points)
    end
    INPUT[] = d
end

nsaved() = length(input()["coords"])
# The first 12 hex digits of the SHA-256 of the input file. Results of different machines are
# only comparable with the same input: a new input has a different Voronoi floe field.
input_checksum() = bytes2hex(open(sha256, INFILE))[1:12]
# Number of copies of the saved floe field, next to each other in x, that n floes need
ntiles(n) = cld(n, nsaved())

const FLOES = Dict{Int, Any}()

"""
    load_floes(n)

`n` floes as a `StructArray{Floe{Float64}}`. More floes than are saved are made by copying
the saved floe field next to itself in x, each copy shifted by `Lx`. Every call returns a new
copy, so steps can change the floes.
"""
function load_floes(n)
    haskey(FLOES, n) && return deepcopy(FLOES[n])
    d = input()
    Lx = d["Lx"]
    floes = Floe{Float64}[]
    for k in 0:n-1
        i = k % nsaved() + 1
        shift = (k ÷ nsaved()) * Lx  # x-offset of this copy of the floe field
        inters = copy(d["interactions"][i])
        inters[:, Int(Subzero.xpoint)] .+= shift
        push!(floes, Floe{Float64}(;
            poly = make_polygon([[[x + shift, y] for (x, y) in d["coords"][i]]]),
            centroid = d["centroid"][i] .+ [shift, 0],
            # relative to the centroid, so they don't need shifting
            x_subfloe_points = copy(d["x_subfloe_points"][i]),
            y_subfloe_points = copy(d["y_subfloe_points"][i]),
            interactions = inters,
            id = k + 1,
            (f => copy(d[string(f)][i]) for f in SAVED_FIELDS)...,
        ))
    end
    FLOES[n] = StructArray(floes)
    return deepcopy(FLOES[n])
end

input_Δt() = Int(input()["Δt"])

"""
    make_model(floes, n)

The model of the input simulation (the docs example `shear_flow.jl` on a larger domain)
with `floes`, for `n` floes: the grid, the periodic domain and the shear ocean flow are
copied next to each other in x like the floes are.
"""
function make_model(floes, n)
    d = input()
    Δgrid, Lx, Ly = d["Δgrid"], d["Lx"], d["Ly"]
    nt = ntiles(n)
    grid = RegRectilinearGrid(; x0 = 0.0, xf = nt * Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
    domain = Domain(;
        north = PeriodicBoundary(North; grid), south = PeriodicBoundary(South; grid),
        east = PeriodicBoundary(East; grid), west = PeriodicBoundary(West; grid),
    )
    # Shear flow in x that changes with y: 0 at the edges of the domain and 0.5 m/s in its
    # middle. The ocean arrays are (Nx + 1) × (Ny + 1), so this is the same for every copy.
    half_ny = fld(grid.Ny, 2)
    u_vec = [range(0, 0.5, length = half_ny); 0.5; range(0.5, 0, length = half_ny)]
    uvels = repeat(transpose(u_vec), outer = (grid.Nx + 1, 1))
    ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)
    atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)
    return Model(; grid, ocean, atmos, domain, floes)
end

# The floe settings of the input simulation. The subfloe point generator is only used for
# new floes (fractures, ridging, welding).
input_floe_settings(grid) = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = input()["npoint_per_cell"]),
)

# The constants of the input simulation. Copying the floe field doesn't change the mean and
# minimum floe area, so this is the same for every n.
function input_constants(floes)
    sqrt_area = sqrt.(floes.area)
    return Constants(; E = 1.5e3 * (mean(sqrt_area) + minimum(sqrt_area)))
end

"""
    ghosted_floes(n)

The floes, model and constants for `n` floes, with ghost floes added as in `timestep_sim!`.
Returns `(floes, model, consts, n_init_floes)`; `floes` aren't `model.floes`.
"""
function ghosted_floes(n)
    floes = load_floes(n)
    consts = input_constants(floes)
    model = make_model(load_floes(0), n)
    Subzero.add_ghosts!(floes, model.domain)
    return floes, model, consts, n
end

"""
    collided_floes(n)

Like `ghosted_floes`, after `timestep_collisions!`: what ridging and rafting starts from.
"""
function collided_floes(n)
    floes, model, consts, n_init = ghosted_floes(n)
    Subzero.timestep_collisions!(floes, n_init, model.domain, consts, input_Δt(),
        CollisionSettings(), Threads.SpinLock())
    return floes, model, consts, n_init
end

"""
    input_simulation(n; kwargs...)

A `Simulation` of the input with `n` floes. `kwargs` are passed on to `Simulation`.
"""
function input_simulation(n; kwargs...)
    floes = load_floes(n)
    model = make_model(floes, n)
    return Simulation(; model, consts = input_constants(floes), Δt = input_Δt(), nΔt = typemax(Int),
        floe_settings = input_floe_settings(model.grid), rng = Xoshiro(1), kwargs...)
end

# Fields of the floes as plain arrays, for the correctness check
floe_fields(floes, fields) = Dict(string(f) => collect(getproperty(floes, f)) for f in fields)
