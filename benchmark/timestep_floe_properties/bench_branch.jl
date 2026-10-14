# Benchmark timestep_floe_properties! on the gpu-version branch (KernelAbstractions kernels).
# Usage: julia -t <nthreads> --project=env-branch bench_branch.jl
# Set BENCH_BACKENDS to choose the backends, e.g. BENCH_BACKENDS=CPU to skip the GPU (for a
# single-thread CPU run), or BENCH_BACKENDS=CPU,AMDGPU. By default, CPU and every GPU
# package (CUDA, AMDGPU) that is installed in the environment and has a working GPU.
include("common.jl")
using Adapt, KernelAbstractions
using Subzero: FixedWidthFloes, timestep_floe_properties!, update_floes!, _log_flags

# GPU packages that can be benchmarked, with their KernelAbstractions backend
const GPU_BACKENDS = Dict("CUDA" => :CUDABackend, "AMDGPU" => :ROCBackend)

const BACKENDS = if haskey(ENV, "BENCH_BACKENDS")
    split(ENV["BENCH_BACKENDS"], ",")
else
    ["CPU"; filter(p -> Base.find_package(p) !== nothing, sort!(collect(keys(GPU_BACKENDS))))]
end
for name in filter(!=("CPU"), BACKENDS)
    haskey(GPU_BACKENDS, name) ||
        error("unknown backend $name, expected CPU, $(join(keys(GPU_BACKENDS), ", "))")
    @eval using $(Symbol(name))
end
gpu_functional(name) = getfield(Main, Symbol(name)).functional()
if haskey(ENV, "BENCH_BACKENDS")
    for name in filter(!=("CPU"), BACKENDS)
        gpu_functional(name) || error("BENCH_BACKENDS contains $name, but $name.functional() is false")
    end
else  # skip GPU packages without a working GPU
    filter!(name -> name == "CPU" || gpu_functional(name), BACKENDS)
end

get_backend(name) = name == "CPU" ? CPU() :
    getfield(getfield(Main, Symbol(name)), GPU_BACKENDS[name])()

function print_device(name)
    if name == "CUDA"
        println(CUDA.name(CUDA.device()), ", CUDA runtime ", CUDA.runtime_version())
    elseif name == "AMDGPU"
        println(AMDGPU.HIP.name(AMDGPU.device()), ", HIP runtime ", AMDGPU.HIP.runtime_version())
    end
end

# Device memory in use, in bytes (by all processes on the device)
function gpu_used_memory(name)
    if name == "CUDA"
        CUDA.total_memory() - CUDA.free_memory()
    elseif name == "AMDGPU"
        free, total = AMDGPU.info()
        total - free
    else
        0
    end
end

# Convert a FixedWidthFloes to Float32 (integer fields are kept as they are)
to_float32(fwf::FixedWidthFloes) = FixedWidthFloes(
    (x -> eltype(x) == Float64 ? Float32.(x) : copy(x)).(getfield.(Ref(fwf), fieldnames(FixedWidthFloes)))...,
)

# Copy all fields of src into dst, to reset dst between samples. This is much faster than
# making a new copy, and doesn't leave garbage that fills up the memory.
function reset!(dst::FixedWidthFloes, src::FixedWidthFloes, backend)
    foreach(f -> copyto!(getfield(dst, f), getfield(src, f)), fieldnames(FixedWidthFloes))
    KernelAbstractions.synchronize(backend)
end

# Free the device memory of floes right away, instead of when they are garbage collected.
# KernelAbstractions.unsafe_free! does nothing for a CuArray (CUDA.jl 6), so use the GPU
# packages' own unsafe_free!.
unsafe_free!(x) = nothing
"CUDA" in BACKENDS && @eval unsafe_free!(x::CUDA.CuArray) = CUDA.unsafe_free!(x)
"AMDGPU" in BACKENDS && @eval unsafe_free!(x::AMDGPU.ROCArray) = AMDGPU.unsafe_free!(x)
free!(floes::FixedWidthFloes) = foreach(f -> unsafe_free!(getfield(floes, f)), fieldnames(FixedWidthFloes))

# What run! does for one call of timestep_floe_properties!, including the conversions and
# copies between host and device
function full_step!(floes, Δt, fs, backend)
    dev_floes = adapt(backend, FixedWidthFloes(floes))
    timestep_floe_properties!(dev_floes, Δt, fs; backend)
    KernelAbstractions.synchronize(backend)
    host_floes = adapt(Array, dev_floes)
    update_floes!(floes, host_floes)
    _log_flags(host_floes.flags, 1, fs)
    free!(dev_floes)
    return
end

function kernels!(dev_floes, Δt, fs, backend)
    timestep_floe_properties!(dev_floes, Δt, fs; backend)
    KernelAbstractions.synchronize(backend)
end

function to_device!(fwf, backend)
    dev = adapt(backend, fwf)
    KernelAbstractions.synchronize(backend)
    free!(dev)
end

foreach(print_device, BACKENDS)
println("gpu-version branch, Julia threads = $(Threads.nthreads())")
fs64, fs32 = floe_settings(Float64), floe_settings(Float32)

# Check that the branch gives the same result as main (saved by bench_main.jl)
if isfile("result_main.jld2")
    ref = load("result_main.jld2", "result")
    nref = load("result_main.jld2", "n")
    for name in BACKENDS
        floes, Δt = load_floes(nref)
        full_step!(floes, Δt, fs64, get_backend(name))
        res = result_fields(floes)
        flat(x) = eltype(x) <: Number ? x : reduce(vcat, vec.(x))
        maxrel = maximum(f -> maximum(abs.(flat(res[f]) .- flat(ref[f]))) /
            max(maximum(abs.(flat(ref[f]))), eps()), keys(res))
        println("$name vs main, n=$nref: all fields isapprox = ",
            all(f -> isapprox(res[f], ref[f]), keys(res)), ", max relative difference = $maxrel")
    end
end

# Every benchmark takes SAMPLES samples (or less if it takes longer than SECONDS). The
# garbage collector runs before every sample (not timed), so that the arrays allocated by
# earlier samples don't fill up the memory.
macro bench(ex, setup = :nothing)
    esc(:(@benchmark($ex, setup = $setup, evals = 1, samples = SAMPLES, seconds = SECONDS,
        gcsample = true)))
end

const GPU_PEAK = Dict(name => 0 for name in BACKENDS)
track_gpu_memory(name) = GPU_PEAK[name] = max(GPU_PEAK[name], gpu_used_memory(name))

for n in SIZES
    floes, Δt = load_floes(n)
    fwf = FixedWidthFloes(floes)
    for name in BACKENDS
        backend = get_backend(name)
        # Kernels only: floes are already on the device. Every sample starts from the
        # floes in src.
        for (FT, fs) in (name == "CPU" ? ((Float64, fs64),) : ((Float64, fs64), (Float32, fs32)))
            src = adapt(backend, FT == Float64 ? fwf : to_float32(fwf))
            d = adapt(backend, deepcopy(src))
            trial = @bench(kernels!($d, $(FT(Δt)), $fs, $backend), reset!($d, $src, $backend))
            track_gpu_memory(name)
            report("branch $name $FT: kernels only", n, trial)
            free!(d)
            FT == Float64 || free!(src)
        end
        # Breakdown of the extra work in run!
        trial = @bench(FixedWidthFloes($floes))
        report("branch $name: FixedWidthFloes(floes)", n, trial)
        trial = @bench(to_device!($fwf, $backend))
        track_gpu_memory(name)
        report("branch $name: adapt to device", n, trial)
        dev = adapt(backend, fwf)
        trial = @bench(adapt(Array, $dev))
        track_gpu_memory(name)
        report("branch $name: adapt to host", n, trial)
        free!(dev)
        # Every sample starts from the floes in fwf. This is the same as resetting them.
        f = deepcopy(floes)
        trial = @bench(update_floes!($f, $fwf))
        report("branch $name: update_floes!", n, trial)
        # Everything together, as in run!
        trial = @bench(full_step!($f, $Δt, $fs64, $backend), update_floes!($f, $fwf))
        track_gpu_memory(name)
        report("branch $name Float64: full step as in run!", n, trial)
    end
end
report_host_memory()
for name in filter(!=("CPU"), BACKENDS)
    @printf("Peak %s memory in use (all processes, measured after each benchmark): %.2f GB\n",
        name, GPU_PEAK[name] / 1e9)
end
