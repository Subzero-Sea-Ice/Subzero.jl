#=
GPU packages to run the tests on. Include this file before loading Subzero or any other
package: adding a GPU package can change the versions of packages that Subzero also uses
(e.g. LLVM), and a package that is already loaded can't be replaced by another version.
=#
import Pkg

# GPU packages that the tests can run on, and the name of their KernelAbstractions backend.
# To support another GPU, add its package here and to `[extras]` and `[compat]` in
# `Project.toml`.
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
#=
A GPU package counts as installed only if it's a dependency of the active (test) environment.
`Base.find_package` isn't enough: it searches the whole `LOAD_PATH`, which by default also
contains the global environment (`@v#.#`), also with TestEnv. A GPU package installed there
would be imported from there instead of being added here, and it would load its own
dependencies (e.g. LLVM) from the global manifest, at versions that may not match the ones
the test environment uses with Subzero. `Pkg.test` isn't affected, because it limits
`LOAD_PATH` to the test environment and the standard library.
=#
let missing_pkgs = filter(pkg -> !haskey(Pkg.project().dependencies, pkg), TEST_GPU_PACKAGES)
    if !isempty(missing_pkgs)
        subzero_dir = dirname(@__DIR__)
        # Don't add test dependencies to Subzero's own Project.toml
        dirname(Base.active_project()) == subzero_dir && error(
            "Run the tests with `Pkg.test()` or TestEnv.jl to install $(join(missing_pkgs, ", ")).",
        )
        # Versions of the packages that are already loaded: these can't change any more
        loaded_versions = Dict(
            id.uuid => pkgversion(mod) for (id, mod) in Base.loaded_modules
            if !isnothing(id.uuid) && !isnothing(pkgversion(mod))
        )
        # The test environment only has `[compat]` entries for the test target, so take
        # the version bounds from Project.toml
        compat = Pkg.TOML.parsefile(joinpath(subzero_dir, "Project.toml"))["compat"]
        Pkg.add([
            Pkg.PackageSpec(; name = pkg, version = Pkg.Types.semver_spec(compat[pkg]))
            for pkg in missing_pkgs
        ])
        #= Adding the GPU packages can change the versions of packages they share with
        Subzero (LLVM, KernelAbstractions, GPUArrays, Adapt, …). If one of those is already
        loaded, the loaded version would be used with the others' new versions. =#
        changed = sort([
            "$(info.name) $(loaded_versions[uuid]) → $(info.version)"
            for (uuid, info) in Pkg.dependencies()
            if haskey(loaded_versions, uuid) && !isnothing(info.version) &&
                info.version != loaded_versions[uuid]
        ])
        isempty(changed) || error(
            "Include gpu_packages.jl before loading other packages: adding ",
            "$(join(missing_pkgs, ", ")) changed the versions of these loaded packages: ",
            join(changed, ", "), ".",
        )
    end
end
for pkg in TEST_GPU_PACKAGES
    @eval import $(Symbol(pkg))
end
