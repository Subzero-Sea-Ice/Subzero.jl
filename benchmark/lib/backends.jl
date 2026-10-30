# KernelAbstractions backends to benchmark the current code on, and helpers for floes on the
# device. Only loaded in the current environment.
#
# Set BENCH_BACKENDS to choose the backends, e.g. BENCH_BACKENDS=CPU to skip the GPU (for a
# single-thread CPU run), or BENCH_BACKENDS=CPU,AMDGPU. By default, CPU and every GPU
# package (CUDA, AMDGPU) that is installed in the environment and has a working GPU.
using Adapt, KernelAbstractions
using Subzero: FixedWidthFloes

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

function device_name(name)
    if name == "CUDA"
        CUDA.name(CUDA.device())
    elseif name == "AMDGPU"
        AMDGPU.HIP.name(AMDGPU.device())
    else
        cpu_name()
    end
end

function print_device(name)
    if name == "CUDA"
        println(device_name(name), ", CUDA runtime ", CUDA.runtime_version())
    elseif name == "AMDGPU"
        println(device_name(name), ", HIP runtime ", AMDGPU.HIP.runtime_version())
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

const GPU_PEAK = Dict(name => 0 for name in BACKENDS)
track_gpu_memory(name) = name == "CPU" || (GPU_PEAK[name] = max(GPU_PEAK[name], gpu_used_memory(name)))
function report_gpu_memory()
    for name in filter(!=("CPU"), BACKENDS)
        @printf("Peak %s memory in use (all processes, measured after each step): %.2f GB\n",
            name, GPU_PEAK[name] / 1e9)
    end
end

# Convert a FixedWidthFloes to Float32 (integer fields are kept as they are)
to_float32(fwf::FixedWidthFloes) = FixedWidthFloes(
    (x -> eltype(x) == Float64 ? Float32.(x) : copy(x)).(getfield.(Ref(fwf), fieldnames(FixedWidthFloes)))...,
)
to_float(::Type{Float64}, fwf::FixedWidthFloes) = fwf
to_float(::Type{Float32}, fwf::FixedWidthFloes) = to_float32(fwf)

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
