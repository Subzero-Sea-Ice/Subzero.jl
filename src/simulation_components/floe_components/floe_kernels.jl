#=
Run per-floe functions on a KernelAbstractions backend (CPU or GPU).
=#

"""
    per_floe_kernel!(backend, workgroup_size)(f, floes, args; ndrange)

KernelAbstractions kernel that calls `f(floes, i, args...)` for each index `i` in `ndrange`.
Launch it with [`launch_per_floe!`](@ref) rather than directly.

## _Positional arguments_
- `f::Function`: per-floe function with the signature `f(floes, i, args...)`. It must be
    GPU-compatible to run on a GPU backend.
- `floes::FixedWidthFloes`: floes on the device of the kernel's backend
- `args::Tuple`: extra arguments passed on to `f`. They must be isbits.

##  _Returns_
- Nothing. `f` updates `floes` in-place.
"""
@kernel function per_floe_kernel!(f, floes, args)
    i = @index(Global)
    f(floes, i, args...)
end

"""
    launch_per_floe!(f, backend, floes, args...; workgroup_size = 512)

Calls `f(floes, i, args...)` for every floe `i` of `floes`, as a kernel on `backend`. `f`
is a plain function, so it can also be called and unit-tested on the CPU without
KernelAbstractions. Functions are isbits, so `f` can be passed to the kernel like any other
argument, and the kernel is compiled once for each `f`.

The kernel runs asynchronously: call `KernelAbstractions.synchronize(backend)` before reading
`floes` on the host. Kernels launched on the same backend run in order.

## _Positional arguments_
- `f::Function`: per-floe function with the signature `f(floes, i, args...)`, for example
    [`calc_strain!`](@ref). It must be GPU-compatible to run on a GPU backend.
- `backend::KernelAbstractions.Backend`: backend to run the kernel on, e.g. `CPU()`
- `floes::FixedWidthFloes`: floes, already moved to `backend` with `adapt`
- `args...`: extra arguments passed on to `f`. They must be isbits.

## _Keyword arguments_
- `workgroup_size::Int`: number of floes per workgroup (Default = 512)

##  _Returns_
- Nothing. `f` updates `floes` in-place.
"""
function launch_per_floe!(f, backend, floes, args...; workgroup_size = 512)
    per_floe_kernel!(backend, workgroup_size)(f, floes, args, ndrange = length(floes.height))
    return
end
