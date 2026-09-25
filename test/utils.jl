#=
GPU packages that the tests can run on, and the name of their KernelAbstractions backend.
To support another GPU, add its package here and to `[extras]` and `[compat]` in
`Project.toml`.
=#
const GPU_BACKENDS = Dict("AMDGPU" => :ROCBackend, "CUDA" => :CUDABackend)

#=
The GPU packages aren't test dependencies, so that not everyone has to install all of them.
Instead, choose the GPUs to test on with a comma-separated list of package names, either as
a test argument, e.g. `Pkg.test(test_args = ["--gpus=AMDGPU"])`, or in the environment
variable `SUBZERO_TEST_GPUS`, e.g. `SUBZERO_TEST_GPUS=AMDGPU`. The test argument takes
precedence. Packages that aren't installed yet are added to the (temporary) test environment.
=#
const TEST_GPU_PACKAGES = let gpus_arg = findlast(startswith("--gpus="), ARGS)
    gpus = isnothing(gpus_arg) ? get(ENV, "SUBZERO_TEST_GPUS", "") :
        chopprefix(ARGS[gpus_arg], "--gpus=")
    filter(!isempty, String.(strip.(split(gpus, ','))))
end
for pkg in TEST_GPU_PACKAGES
    haskey(GPU_BACKENDS, pkg) || error(
        "Unknown GPU package $pkg, choose from $(join(keys(GPU_BACKENDS), ", ")).",
    )
end
let missing_pkgs = filter(pkg -> isnothing(Base.find_package(pkg)), TEST_GPU_PACKAGES)
    if !isempty(missing_pkgs)
        # Don't add test dependencies to Subzero's own Project.toml
        dirname(Base.active_project()) == pkgdir(Subzero) && error(
            "Run the tests with `Pkg.test()` or TestEnv.jl to install $(join(missing_pkgs, ", ")).",
        )
        import Pkg
        # The test environment only has `[compat]` entries for the test target, so take
        # the version bounds from Project.toml
        compat = Pkg.TOML.parsefile(joinpath(pkgdir(Subzero), "Project.toml"))["compat"]
        Pkg.add([
            Pkg.PackageSpec(; name = pkg, version = Pkg.Types.semver_spec(compat[pkg]))
            for pkg in missing_pkgs
        ])
    end
end
for pkg in TEST_GPU_PACKAGES
    @eval import $(Symbol(pkg))
end

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

# JLArrays (0.3.3) doesn't define `synchronize` for its backend. Its kernels run
# synchronously, so there is nothing to wait for.
Subzero.KernelAbstractions.synchronize(::JLArrays.JLBackend) = nothing

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
