#=
Compare a whole simulation on every backend with the results of the original CPU
implementation (Subzero at ae96ff2, before the GPU port), which
inputs/make_simulation_reference.jl saved in inputs/simulation_reference.jld2. The
simulation (in inputs/simulation_reference_setup.jl) is a small version of the shear flow
example, with periodic boundaries, collisions, and coupling, and starts from the floes saved
with the reference results.

Round-off differences grow over time through the collisions, so later timesteps are compared
more loosely. The tolerances are about 100 times the largest relative difference of a floe's
value that was measured on the CPU, JLArrays, and AMDGPU backends when the reference was made
(1e-15, 2e-11, 5e-10, and 1e-7 after timesteps 1, 10, 50, and 100).
=#
include("inputs/simulation_reference_setup.jl")

const REFERENCE_RTOL = Dict(1 => 1e-13, 10 => 1e-9, 50 => 1e-7, 100 => 1e-5)

@testset "Simulation matches the original implementation on $backend" for backend in test_backends()
    saved = load("inputs/simulation_reference.jld2")
    sim = reference_simulation(reference_floes_from_dict(saved["floes"]); backend)
    results = run_reference_simulation(sim)
    @testset "after timestep $tstep" for tstep in REFERENCE_TIMESTEPS
        got, expected = results[tstep], saved["results"][string(tstep)]
        @test got["id"] == expected["id"]
        if got["id"] == expected["id"]
            # the fields that differ, with the floes in which they differ
            rtol = REFERENCE_RTOL[tstep]
            differing = filter(!isempty ∘ last, [
                f => _not_approx(got[f], expected[f]; rtol)
                for f in [string.(filter(!=(:id), REFERENCE_RESULT_FIELDS))..., "coords"]
            ])
            @test isempty(differing)
        end
    end
end
