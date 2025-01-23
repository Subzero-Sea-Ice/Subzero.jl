using Test, LinearAlgebra
import GeometryOps.GeoInterface as GI
# Create test floe
poly = GI.Polygon([[(0.0, 0.0), (0.0, 100.0), (100.0, 0.0), (0.0, 0.0)]])
hmean, Δh = 0.5, 0.0
stress_accum = [0.3 0.5; 0.7 0.9]
floe_settings = FloeSettings(min_floe_area = 2.5e3)
floe = Floe(poly, hmean, Δh; stress_accum, floe_settings)

# # Test DecayAreaScaledCalculator
c1 = DecayAreaScaledCalculator()
c2 = (@test_logs (:warn, Subzero.DECAY_ARG_WARNING) DecayAreaScaledCalculator(λ = 1.1, α = 2.0))
c3 = DecayAreaScaledCalculator(Float32, λ = 0.2, α = 2)

# Test constructors
@test c1.λ == 0.2 && c1.α == 0
@test c1 isa DecayAreaScaledCalculator{Float64}
@test c2.λ == 0.2 && c2.α == 2.0
@test c3.λ == 0.2f0 && c3.α == 2.0f0
@test c3 isa DecayAreaScaledCalculator{Float32}

# Test _update_stress_accum!
curr_stress = [0.2 0.4; 0.6 0.8]
new_stress_accum = 0.8 * stress_accum .+ 0.2 * curr_stress
Subzero._update_stress_accum!(c1, curr_stress, floe)
@test all(floe.stress_accum .== new_stress_accum)

# Test _scale_principal_stress!
σvals = eigvals(floe.stress_accum)
new_σvals = 4σvals
Subzero._scale_principal_stress!(c2, σvals, floe, floe_settings)
@test all(σvals .== new_σvals)


# # Test DamageStressCalculator
# TODO: Add tests for DamageStressCalculator when it is implemented
@test_throws ErrorException DamageStressCalculator()
