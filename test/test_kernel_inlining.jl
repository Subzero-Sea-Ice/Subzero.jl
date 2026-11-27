#=
The kernels must not copy the floes into each thread's memory. If a function that takes the
floes, or one of their arrays, isn't inlined into the kernel, the GPU compiler passes them to
it by reference, and so first copies the whole `FixedWidthFloes` (over 1 kB) into each thread's
local memory. That made the kernels several times slower on CUDA.

The test compiles the kernels for the CPU backend (POCL), which uses the same Julia-level
inlining as the GPU backends, and looks at their optimized Julia IR. So it runs without a GPU.
It doesn't see what LLVM does afterwards: on CUDA, the exact check is that the PTX has no
`__local_depot` (`CUDA.@device_code_ptx`).
=#

# The type of argument `a` of a statement in the optimized code `ci`
_argtype(a::Core.SSAValue, ci) = Core.Compiler.widenconst(ci.ssavaluetypes[a.id])
_argtype(a::Core.Argument, ci) = Core.Compiler.widenconst(ci.slottypes[a.n])
_argtype(a::GlobalRef, ci) = typeof(getglobal(a.mod, a.name))
_argtype(a::QuoteNode, ci) = typeof(a.value)
_argtype(a, ci) = typeof(a)

# Whether a value of type T is or holds an array, e.g. the floes or a polygon view of them
_holds_array(T) = T <: AbstractArray ||
    (isconcretetype(T) && isstructtype(T) && any(_holds_array, fieldtypes(T)))

#=
The calls that weren't inlined into the kernels that `f()` launches on the CPU backend and that
take an array, as a Dict from the kernel (for `per_floe_kernel!`, the per-floe function) to the
names of the called functions. Functions whose names start with `throw` are left out: they are
never inlined, and on the GPU they don't make the kernel copy its arguments.
=#
function _not_inlined_calls(f)
    kernels = Subzero.KernelAbstractions.POCL.@device_code_typed f()
    result = Dict{String, Vector{String}}()
    for (job, codes) in kernels, (ci, _) in codes
        calls = String[]
        for stmt in ci.code
            Meta.isexpr(stmt, :invoke) || continue
            mi = stmt.args[1]
            mi isa Core.CodeInstance && (mi = mi.def)
            name = string(mi.def.name)
            startswith(name, "throw") && continue
            any(a -> _holds_array(_argtype(a, ci)), stmt.args[3:end]) && push!(calls, name)
        end
        isempty(calls) && continue
        types = job.source.specTypes.parameters
        kernel = string(nameof(types[1].instance))
        # per_floe_kernel!(ctx, f, floes, args): name it after f
        if occursin("per_floe_kernel", kernel)
            kernel = string(nameof(types[3].instance))
        end
        result[kernel] = unique!(calls)
    end
    return result
end

@testset "Kernels inline the functions that take the floes" begin
    backend = Subzero.KernelAbstractions.CPU()
    #= TODO: timestep_floe! doesn't inline calc_stress!, calc_strain!, limit_collision_force!,
    thermodynamic_growth!, update_ice_coordinates! and update_velocities!, so on CUDA each
    thread copies the floes into local memory. To be fixed in a separate PR. =#
    @testset "timestep_floe_properties! with $FT" for FT in (Float64, Float32)
        floe_settings = FloeSettings(FT)
        floes = Subzero.FixedWidthFloes(_make_timestep_test_floes(FT, floe_settings))
        @test_broken isempty(_not_inlined_calls(
            () -> Subzero.timestep_floe_properties!(floes, 10, floe_settings; backend)))
    end
    @testset "add_ghosts! and remove_ghosts! with $FT" for FT in (Float64, Float32)
        grid = RegRectilinearGrid(FT; x0 = 0.0, xf = 1e4, y0 = 0.0, yf = 1e4, Δx = 1e3, Δy = 1e3)
        domain = Domain(; north = PeriodicBoundary(North, FT; grid),
            south = PeriodicBoundary(South, FT; grid), east = PeriodicBoundary(East, FT; grid),
            west = PeriodicBoundary(West, FT; grid))
        floes = _ghost_test_floes(FT)
        # with room for the ghosts, so that every kernel of both passes runs
        fwf = Subzero.FixedWidthFloes(floes; capacity = 2length(floes))
        @test isempty(_not_inlined_calls(function ()
            ghosted, _ = Subzero.add_ghosts!(fwf, domain, backend)
            Subzero.remove_ghosts!(ghosted, length(floes), backend)
        end))
    end
end
