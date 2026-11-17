#=
Run an example with one version of Subzero, from the saved floes, and save what is needed to
draw it every `every` timesteps:
    julia --project=env-<version> run_example.jl <example> <version> <nsteps> <every> <floes.jld2> <out.jld2>
`version` is `original` or `current`. The current version runs on a GPU if the environment
has a GPU package (AMDGPU or CUDA) that works, and else on the CPU backend. Set
`VIS_BACKEND=CPU` to run it on the CPU backend anyway.

For each frame, before the timestep (as the documentation's FloeOutputWriter):
- `polys`: the exterior rings of the floes and of their periodic ghosts, so that floes that
  cross a periodic boundary are drawn on both sides, as in the documentation's animations;
- `ids`, `centroids`: the ids and centroids of the floes without ghosts, to compare the
  versions floe by floe;
- `walls`: the positions (`val`) of the north, south, east, and west boundaries, which
  change for moving boundaries.
Everything is saved as plain arrays, so that the other version of Subzero can read it.
=#
using Subzero, JLD2, Logging, Random, Statistics, StructArrays
import Pkg
const GI = Subzero.GI
include(joinpath(@__DIR__, "examples.jl"))

name, version = ARGS[1], ARGS[2]
nsteps, every = parse(Int, ARGS[3]), parse(Int, ARGS[4])
floes_file, out_file = ARGS[5], ARGS[6]

# The first GPU package of the environment that works, else the CPU backend
function current_backend()
    get(ENV, "VIS_BACKEND", "") == "CPU" && return Subzero.KernelAbstractions.CPU()
    for (pkg, backend) in (("AMDGPU", :ROCBackend), ("CUDA", :CUDABackend))
        haskey(Pkg.project().dependencies, pkg) || continue
        mod = Base.require(Main, Symbol(pkg))
        Base.invokelatest(mod.functional) || continue
        return Base.invokelatest(getfield(mod, backend))
    end
    return Subzero.KernelAbstractions.CPU()
end

kwargs = version == "current" ? (; backend = current_backend()) : (;)
backend_name = version == "current" ? string(nameof(typeof(kwargs.backend))) : "CPU (threads)"
floes = floes_from_dict(load(floes_file, "floes"))
sim = example_simulation(name, floes, nsteps; kwargs...)

frames = (tsteps = Int[], polys = Vector{Vector{NTuple{2, Float64}}}[],
    ids = Vector{Int}[], centroids = Vector{NTuple{2, Float64}}[],
    walls = NTuple{4, Float64}[])
function record_frame!(frames, sim, tstep)
    floes = sim.model.floes
    with_ghosts = deepcopy(floes)
    Subzero.add_ghosts!(with_ghosts, sim.model.domain)
    push!(frames.tsteps, tstep)
    push!(frames.polys, ring_coords(with_ghosts.poly))
    push!(frames.ids, collect(floes.id))
    push!(frames.centroids, [(Float64(c[1]), Float64(c[2])) for c in floes.centroid])
    d = sim.model.domain
    push!(frames.walls, Float64.((d.north.val, d.south.val, d.east.val, d.west.val)))
    return
end

elapsed = @elapsed with_logger(NullLogger()) do
    for tstep in 0:nsteps
        tstep % every == 0 && record_frame!(frames, sim, tstep)
        tstep < nsteps && timestep_sim!(sim, tstep)
    end
end
jldsave(out_file; frames.tsteps, frames.polys, frames.ids, frames.centroids, frames.walls,
    backend = backend_name, commit = get(ENV, "VIS_COMMIT", "unknown"),
    threads = Threads.nthreads(), elapsed)
println("$name, $version on $backend_name: $nsteps timesteps in $(round(elapsed, digits = 1)) s, ",
    "$(length(sim.model.floes)) floes at the end")
