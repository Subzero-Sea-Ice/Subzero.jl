#=
Helpers for the tests of the GeometryOpsGPU submodule (`test_geometry_ops_gpu*.jl`). They use
only GeometryOpsGPU, GeometryOps and the packages they build on, not Subzero, so that the
tests can move with GeometryOpsGPU if it becomes a package of its own. The one exception is
the list of backends to run kernels on, `test_backends()` in `utils.jl`.
=#
import Adapt, CoordinateTransformations, Extents, Rotations
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

# A global constant, so that kernels defined in testsets can use it without capturing it
const GOG = Subzero.GeometryOpsGPU

#= Polygons for testing: random polygons around the origin, so that most pairs overlap.
Star-shaped polygons with 3 to 9 vertices, and rectangles and L-shapes on a 250 m grid, which
share edges and vertices with each other. Some are clockwise. =#
function _geom_test_polys(::Type{FT}; n = 150, rng = Xoshiro(1)) where FT
    polys = Vector{Vector{Tuple{FT, FT}}}()
    while length(polys) < n
        if rand(rng) < 0.4
            nv = rand(rng, 3:9)
            θs = sort(2π .* rand(rng, nv))
            rs = rand(rng) < 0.5 ? fill(1e3, nv) : 1e3 .* (0.3 .+ rand(rng, nv))
            c = 500 .* (rand(rng, 2) .- 0.5)
            push!(polys, [(FT(c[1] + r * cos(θ)), FT(c[2] + r * sin(θ))) for (r, θ) in zip(rs, θs)])
        else
            g = FT.(250 .* sort(rand(rng, -4:4, 3)))
            h = FT.(250 .* sort(rand(rng, -4:4, 3)))
            (allunique(g) && allunique(h)) || continue
            push!(polys, rand(rng) < 0.5 ?
                [(g[1], h[1]), (g[3], h[1]), (g[3], h[3]), (g[1], h[3])] :
                [(g[1], h[1]), (g[3], h[1]), (g[3], h[2]), (g[2], h[2]), (g[2], h[3]), (g[1], h[3])])
        end
        rand(rng) < 0.3 && reverse!(polys[end])
    end
    return polys
end
#= The polygons as `GOG.ArrayPolygon`s, closed rings in rows of a 3D array, padded with
zeros. Returns the array and the polygons. =#
function _array_polys(polys::Vector{Vector{Tuple{FT, FT}}}) where FT
    coords = zeros(FT, length(polys), maximum(length, polys) + 1, 2)
    apolys = map(enumerate(polys)) do (i, p)
        for (j, (x, y)) in enumerate([p; [p[1]]])
            coords[i, j, 1], coords[i, j, 2] = x, y
        end
        GOG.ArrayPolygon(coords, i, length(p) + 1)
    end
    return coords, apolys
end
_gi_poly(p) = GI.Polygon([GI.LinearRing([p; [p[1]]])])

#= The exterior rings of the floe shapes in test/inputs/floe_shapes.jld2 (7 to 591 vertices,
coordinates up to ±65 km), as vectors of (x, y), without repeating the first point =#
_floe_shape_rings() = jldopen(joinpath(@__DIR__, "inputs", "floe_shapes.jld2")) do f
    [[(c[1], c[2]) for c in p[1]] for p in vec(f["floe_vertices"])]
end

# A test point for polygon p: a vertex, a point along an edge, a grid point, or a random point
function _geom_test_points(rng, p::Vector{Tuple{FT, FT}}, polys) where FT
    u = rand(rng)
    pt = if u < 0.3
        rand(rng, rand(rng, polys))
    elseif u < 0.5
        (FT(250 * rand(rng, -4:4)), FT(125 * rand(rng, -10:10)))
    elseif u < 0.7
        v, w = rand(rng, p), rand(rng, p)
        v .+ FT(rand(rng, 0:4) / 4) .* (w .- v)
    else
        (FT(3e3 * rand(rng) - 1.5e3), FT(3e3 * rand(rng) - 1.5e3))
    end
    return pt
end

#= The name of the kernel of a compilation job. A kernel that runs a function given as its
first argument, such as `_gog_test_kernel!(f, args)`, is named after that function. =#
function _kernel_name(job)
    # the kernel, the KernelAbstractions context, and the kernel's arguments
    kernel, _, args... = job.source.specTypes.parameters
    f = isempty(args) ? kernel : args[1]
    return string(nameof(f <: Function && isdefined(f, :instance) ? f.instance : kernel.instance))
end

#=
The kernels that `f()` launches on the CPU backend whose final LLVM IR has a `double` or an
`fpext`, i.e. that compute in Float64. Unlike counting the Float64 values in the typed code
of one function, this sees everything the kernel calls, inlined or not.
=#
function _float64_kernels(f)
    jobs = keys(KA.POCL.@device_code_typed f())
    return unique!([_kernel_name(job) for job in jobs if _computes_in_float64(job)])
end

# Whether the final LLVM IR of a compilation job has a `double` or an `fpext`
_computes_in_float64(job) = occursin(r"\bdouble\b|\bfpext\b",
    sprint(io -> KA.POCL.GPUCompiler.code_llvm(io, job; dump_module = true)))

#= Whether the first kernel that `launch()` launches on the CPU backend computes in Float64,
as `_float64_kernels` checks, but without compiling the kernel to machine code or running it:
generating its LLVM IR takes seconds, where POCL can take minutes to compile a large kernel.
The kernel must not have been compiled already. =#
struct _StopCompile <: Exception end
function _kernel_computes_in_float64(launch)
    job = Ref{Any}()
    try
        KA.POCL.GPUCompiler.with_compile_hook(launch, j -> (job[] = j; throw(_StopCompile())))
    catch e
        e isa _StopCompile || rethrow()
    end
    isassigned(job) || error("no kernel was compiled: it may have been compiled already")
    return _computes_in_float64(job[])
end

#= Run `launch()`. With Float32 on the CPU backend, also check that the kernels it launches
don't compute in Float64, so that they can run on GPUs without Float64, such as Apple GPUs
(Metal) and some Intel GPUs. =#
function _launch_float32_checked(launch, backend, FT)
    if FT == Float32 && backend isa typeof(KA.CPU())
        @test isempty(_float64_kernels(launch))
    else
        launch()
    end
end

# The points of a GOG.ArrayPolygon, and whether its padding after them is zero
_gog_points(p) = [GOG.getpoint(p, j) for j in 1:GOG.npoint(p)]
_gog_padding_is_zero(p) = all(iszero, p.coords[p.i, (GOG.npoint(p) + 1):end, :])

# The rectangle `ext` as a GeoInterface polygon, for GeometryOps
_gi_rect(ext) = (((x1, x2), (y1, y2)) = (ext.X, ext.Y);
    GI.Polygon([GI.LinearRing([(x1, y1), (x1, y2), (x2, y2), (x2, y1), (x1, y1)])]))

#= Run `f(i, args...)` for `i in 1:n` in a kernel on `backend`, and wait for it.
Arguments must be on the device and isbits, e.g. arrays moved with
`Adapt.adapt(backend, x)`. Let KernelAbstractions choose the workgroup size because
a fixed one can exceed the limit of a large kernel (in terms of estimated stack usage).
This depends on the machine (e.g. PoCL on github CI only allowed 224 work-items for a
large kernel, resulting in an "ArgumentError: Work-group size has more than 224
work-items, the maximum for this kernel" error). =#
@kernel function _gog_test_kernel!(f, args)
    i = @index(Global)
    @inline f(i, args...)
end
function _gog_launch!(f, backend, n, args...)
    _gog_test_kernel!(backend)(f, args; ndrange = n)
    KA.synchronize(backend)
    return
end
