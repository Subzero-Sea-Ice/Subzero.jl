# Benchmark timestep_floe_properties! on the gpu-version branch (KernelAbstractions kernels).
# Usage: julia -t <nthreads> --project=env-branch bench_branch.jl
# Set BENCH_BACKENDS=CPU to skip CUDA (e.g. for a single-thread CPU run).
include("common.jl")
using CUDA, Adapt, KernelAbstractions
using Subzero: FixedWidthFloes, timestep_floe_properties!, update_floes!, _log_flags

const BACKENDS = split(get(ENV, "BENCH_BACKENDS", "CPU,CUDA"), ",")

# Convert a FixedWidthFloes to Float32 (integer fields are kept as they are)
to_float32(fwf::FixedWidthFloes) = FixedWidthFloes(
    (x -> eltype(x) == Float64 ? Float32.(x) : copy(x)).(getfield.(Ref(fwf), fieldnames(FixedWidthFloes)))...,
)

# What run! does for one call of timestep_floe_properties!, including the conversions and
# copies between host and device
function full_step!(floes, Δt, fs, backend)
    dev_floes = adapt(backend, FixedWidthFloes(floes))
    timestep_floe_properties!(dev_floes, Δt, fs; backend)
    KernelAbstractions.synchronize(backend)
    host_floes = adapt(Array, dev_floes)
    update_floes!(floes, host_floes)
    _log_flags(host_floes.flags, 1, fs)
    return
end

function kernels!(dev_floes, Δt, fs, backend)
    timestep_floe_properties!(dev_floes, Δt, fs; backend)
    KernelAbstractions.synchronize(backend)
end

CUDA.functional() && println(CUDA.name(CUDA.device()), ", CUDA runtime ", CUDA.runtime_version())
println("gpu-version branch, Julia threads = $(Threads.nthreads())")
fs64, fs32 = floe_settings(Float64), floe_settings(Float32)

# Check that the branch gives the same result as main (saved by bench_main.jl)
if isfile("result_main.jld2")
    ref = load("result_main.jld2", "result")
    nref = load("result_main.jld2", "n")
    for name in BACKENDS
        floes, Δt = load_floes(nref)
        full_step!(floes, Δt, fs64, name == "CUDA" ? CUDABackend() : CPU())
        res = result_fields(floes)
        flat(x) = eltype(x) <: Number ? x : reduce(vcat, vec.(x))
        maxrel = maximum(f -> maximum(abs.(flat(res[f]) .- flat(ref[f]))) /
            max(maximum(abs.(flat(ref[f]))), eps()), keys(res))
        println("$name vs main, n=$nref: all fields isapprox = ",
            all(f -> isapprox(res[f], ref[f]), keys(res)), ", max relative difference = $maxrel")
    end
end
for n in SIZES
    floes, Δt = load_floes(n)
    fwf = FixedWidthFloes(floes)
    for name in BACKENDS
        backend = name == "CUDA" ? CUDABackend() : CPU()
        # Kernels only: floes are already on the device
        trial = @benchmark(kernels!(d, $Δt, $fs64, $backend),
            setup = (d = adapt($backend, deepcopy($fwf))), evals = 1, seconds = SECONDS)
        report("branch $name Float64: kernels only", n, trial)
        if name == "CUDA"
            fwf32 = to_float32(fwf)
            trial = @benchmark(kernels!(d, $Δt, $fs32, $backend),
                setup = (d = adapt($backend, deepcopy($fwf32))), evals = 1, seconds = SECONDS)
            report("branch $name Float32: kernels only", n, trial)
        end
        # Breakdown of the extra work in run!
        trial = @benchmark(FixedWidthFloes($floes), seconds = SECONDS)
        report("branch $name: FixedWidthFloes(floes)", n, trial)
        trial = @benchmark((adapt($backend, $fwf); KernelAbstractions.synchronize($backend)), seconds = SECONDS)
        report("branch $name: adapt to device", n, trial)
        dev = adapt(backend, fwf)
        trial = @benchmark(adapt(Array, $dev), seconds = SECONDS)
        report("branch $name: adapt to host", n, trial)
        trial = @benchmark(update_floes!(f, $fwf),
            setup = (f = deepcopy($floes)), evals = 1, seconds = SECONDS)
        report("branch $name: update_floes!", n, trial)
        # Everything together, as in run!
        trial = @benchmark(full_step!(f, $Δt, $fs64, $backend),
            setup = (f = deepcopy($floes)), evals = 1, seconds = SECONDS)
        report("branch $name Float64: full step as in run!", n, trial)
    end
end
