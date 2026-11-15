# The GPU packages, if gpu_packages.jl wasn't included before loading Subzero
isdefined(@__MODULE__, :TEST_GPU_PACKAGES) || include("gpu_packages.jl")

#=
KernelAbstractions backends to run the kernel tests on: the CPU, JLArrays' reference GPU
backend (which runs on the CPU, so the GPU code path is tested without a GPU), and the GPUs
chosen with `SUBZERO_TEST_GPUS`.
=#
const TEST_BACKENDS = let backends = Any[Subzero.KernelAbstractions.CPU(), JLArrays.JLBackend()]
    for pkg in TEST_GPU_PACKAGES
        gpu = getfield(@__MODULE__, Symbol(pkg))
        Base.invokelatest(gpu.functional) || error(
            "$pkg can't use a GPU on this computer, see $pkg.jl's documentation.",
        )
        push!(backends, Base.invokelatest(getfield(gpu, GPU_BACKENDS[pkg])))
    end
    backends
end
test_backends() = TEST_BACKENDS

# JLArrays 0.3 doesn't define `synchronize` for its backend (0.5 does). Its kernels run
# synchronously, so there is nothing to wait for.
if !hasmethod(Subzero.KernelAbstractions.synchronize, Tuple{JLArrays.JLBackend})
    Subzero.KernelAbstractions.synchronize(::JLArrays.JLBackend) = nothing
end

#=
Helpers to compare the values of many floes in one test: a failure then lists the floes (or
cases) that differ, instead of being one of many tests per floe.

`_not_approx(a, b; kwargs...)` gives the indices of the elements of `a` and `b` that aren't
`isapprox`. Each element is compared as a whole, as `isapprox(a[i], b[i]; kwargs...)`, e.g.
a floe's centroid or stress, and elements of different sizes aren't approx. This is stricter
than `isapprox(a, b)`, which compares the norm of the difference of the whole arrays. Gives
`[:size]` if `a` and `b` have different sizes.
`_failing(f, is)` gives the elements `i` of `is` for which `f(i)` is false.
=#
function _not_approx(a, b; kwargs...)
    size(a) == size(b) || return [:size]
    return findall(map((x, y) -> !(size(x) == size(y) && isapprox(x, y; kwargs...)), a, b))
end
_failing(f, is) = filter(!f, collect(is))

# The vertices of floe i of a FixedWidthFloes, and whether the padding after them is zero
_fwf_points(fwf, i) = [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:fwf.n_points[i]]
_fwf_padding_is_zero(fwf, i) = all(iszero, fwf.poly[i, (fwf.n_points[i] + 1):end, :])

#=
Run the per-floe function `f` on a copy of `fwf::FixedWidthFloes` on `backend` with
`launch_per_floe!` and return the result as a `FixedWidthFloes` on the CPU.
=#
function _launch_on(f, backend, fwf, args...)
    dev_floes = Subzero.adapt(backend, fwf)
    Subzero.launch_per_floe!(f, backend, dev_floes, args...)
    Subzero.KernelAbstractions.synchronize(backend)
    return Subzero.adapt(Array, dev_floes)
end

#=
    find_poly_coords(poly)

Syntactic sugar for to find a polygon's coordinates
Input:
    poly    <Polygon>
Output:
    <PolyVec> representing the floe's coordinates xy plane
=#
find_poly_coords(poly) = GI.coordinates(poly)

# Make a copy of given coordinates and translate by given deltas, returning new coordiantes
function translate_coords(coords, Δx, Δy)
    new_coords = [[Vector{Float64}(undef, 2) for _ in eachindex(coords[1])]]
    for i in eachindex(coords[1])
        new_coords[1][i][1] = coords[1][i][1] + Δx
        new_coords[1][i][2] = coords[1][i][2] + Δy 
    end
    return new_coords
end

# Translate each of the given coodinates by given deltas in place
function translate_coords!(coords, Δx, Δy)
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            coords[i][j][1] += Δx
            coords[i][j][2] += Δy
        end
    end
    return
end

#=
Deterministic set of floes for testing `timestep_floe_properties!`. Floes have different
numbers of vertices and non-zero velocities, forces, and interactions. Some floes trigger
the height limit, collision force reduction, velocity adjustment, and ξ limit.
=#
function _make_timestep_test_floes(::Type{FT}, floe_settings; n = 10) where FT
    rng = Xoshiro(42)
    floes = StructArray{Floe{FT}}(undef, 0)
    for k in 1:n
        # convex polygon with k + 3 vertices around a random center
        nverts = k + 3
        cx, cy = 1e4 .* rand(rng, 2)
        r = 500 + 1500rand(rng)
        θs = sort(2π .* rand(rng, nverts))
        ring = [[cx + r * cos(θ), cy + r * sin(θ)] for θ in θs]
        push!(ring, ring[1])
        height = 0.5 + rand(rng)
        floe = Floe{FT}([ring], height; floe_settings, rng)
        floe.u, floe.v = 0.2 .* (rand(rng, 2) .- 0.5)
        floe.ξ = 1e-6 * (rand(rng) - 0.5)
        floe.p_dxdt, floe.p_dydt = 0.2 .* (rand(rng, 2) .- 0.5)
        floe.p_dαdt = 1e-6 * (rand(rng) - 0.5)
        floe.p_dudt, floe.p_dvdt = 1e-4 .* (rand(rng, 2) .- 0.5)
        floe.p_dξdt = 1e-9 * (rand(rng) - 0.5)
        floe.fxOA, floe.fyOA = 1e6 .* (rand(rng, 2) .- 0.5)
        floe.trqOA = 1e8 * (rand(rng) - 0.5)
        floe.hflx_factor = 1e-3 * (rand(rng) - 0.5)
        floe.collision_force = 1e6 .* (rand(rng, 1, 2) .- 0.5)
        floe.collision_trq = 1e8 * (rand(rng) - 0.5)
        floe.stress_accum = 1e3 .* (rand(rng, 2, 2) .- 0.5)
        floe.stress_instant = 1e3 .* (rand(rng, 2, 2) .- 0.5)
        floe.strain = 1e-6 .* (rand(rng, 2, 2) .- 0.5)
        floe.damage = rand(rng)
        if isodd(k)  # interactions, with unused (padding) rows for some floes
            ninters = k ÷ 2 + 1
            inters = zeros(FT, ninters + (k % 3), 7)
            inters .= 1e3 .* rand(rng, size(inters)...)
            inters[:, xpoint] .+= cx
            inters[:, ypoint] .+= cy
            floe.interactions = inters
            floe.num_inters = ninters
        end
        push!(floes, floe)
    end
    # edge cases
    floes.height[1] = 2floe_settings.max_floe_height  # height limit
    floes.collision_force[2] .= [1e3 -1e3] .* floes.mass[2]  # force reduction
    floes.fxOA[3] = 10 * floes.mass[3]  # velocity adjustment (frac != 1)
    floes.fyOA[4] = -10 * floes.mass[4]
    floes.fxOA[4] = 10 * floes.mass[4]
    floes.trqOA[5] = 1e3 * floes.moment[5]  # ξ limit
    return floes
end

# The name of the kernel of a compilation job, or for `per_floe_kernel!`, the per-floe function
function _kernel_name(job)
    types = job.source.specTypes.parameters
    kernel = string(nameof(types[1].instance))
    # per_floe_kernel!(ctx, f, floes, args): name it after f
    if occursin("per_floe_kernel", kernel)
        kernel = string(nameof(types[3].instance))
    end
    return kernel
end

#=
The kernels that `f()` launches on the CPU backend whose final LLVM IR has a `double` or an
`fpext`, i.e. that compute in Float64. Unlike counting the Float64 values in the typed code
of one function, this sees everything the kernel calls, inlined or not.
=#
function _float64_kernels(f)
    jobs = keys(Subzero.KernelAbstractions.POCL.@device_code_typed f())
    llvm(job) = sprint(io -> Subzero.KernelAbstractions.POCL.GPUCompiler.code_llvm(io, job;
        dump_module = true))
    return unique!([_kernel_name(job) for job in jobs if
        occursin(r"\bdouble\b|\bfpext\b", llvm(job))])
end

#= Polygons for testing the device geometry functions: random polygons around the origin,
so that most pairs overlap. Star-shaped polygons with 3 to 9 vertices, and rectangles and
L-shapes on a 250 m grid, which share edges and vertices with each other. Some are
clockwise. =#
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
# The polygons as `GOG.ArrayPolygon`s: closed rings in a 3D array, with the layout of
# FixedWidthFloes.poly
function _array_polys(polys::Vector{Vector{Tuple{FT, FT}}}) where FT
    coords = zeros(FT, length(polys), maximum(length, polys) + 1, 2)
    apolys = map(enumerate(polys)) do (i, p)
        for (j, (x, y)) in enumerate([p; [p[1]]])
            coords[i, j, 1], coords[i, j, 2] = x, y
        end
        Subzero.GeometryOpsGPU.ArrayPolygon(coords, i, length(p) + 1)
    end
    return coords, apolys
end
_gi_poly(p) = GI.Polygon([GI.LinearRing([p; [p[1]]])])

#= The exterior rings of the floe shapes in test/inputs/floe_shapes.jld2 (7 to 591 vertices,
coordinates up to ±65 km), as vectors of (x, y), without repeating the first point =#
_floe_shape_rings() = jldopen(joinpath(@__DIR__, "inputs", "floe_shapes.jld2")) do f
    [[(c[1], c[2]) for c in p[1]] for p in vec(f["floe_vertices"])]
end
