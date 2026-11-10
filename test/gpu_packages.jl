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
let missing_pkgs = filter(pkg -> isnothing(Base.find_package(pkg)), TEST_GPU_PACKAGES)
    if !isempty(missing_pkgs)
        subzero_dir = dirname(@__DIR__)
        # Don't add test dependencies to Subzero's own Project.toml
        dirname(Base.active_project()) == subzero_dir && error(
            "Run the tests with `Pkg.test()` or TestEnv.jl to install $(join(missing_pkgs, ", ")).",
        )
        loaded = filter(id -> id.name in ("Subzero", "LLVM"), keys(Base.loaded_modules))
        isempty(loaded) || error(
            "Include gpu_packages.jl before loading $(join(sort([id.name for id in loaded]), " and ")): ",
            "adding $(join(missing_pkgs, ", ")) may need other versions of the packages they use.",
        )
        # The test environment only has `[compat]` entries for the test target, so take
        # the version bounds from Project.toml
        compat = Pkg.TOML.parsefile(joinpath(subzero_dir, "Project.toml"))["compat"]
        Pkg.add([
            Pkg.PackageSpec(; name = pkg, version = Pkg.Types.semver_spec(compat[pkg]))
            for pkg in missing_pkgs
        ])
    end
end
for pkg in TEST_GPU_PACKAGES
    @eval import $(Symbol(pkg))
end
