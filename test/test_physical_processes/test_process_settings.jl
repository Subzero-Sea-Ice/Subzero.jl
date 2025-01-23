@testset "ProcessInfo" begin
    @testset "FloeSettings" begin
        # Test default
        default_info = FloeSettings()
        @test default_info.ρi == 920 &&
            default_info.min_floe_area == 1e6 &&
            default_info.min_floe_height == 0.1 &&
            default_info.max_floe_height == 10 &&
            default_info.min_aspect_ratio == 0.05 &&
            default_info.stress_calculator isa DecayAreaScaledCalculator &&
            default_info.subfloe_point_generator isa MonteCarloPointsGenerator
        @test default_info.ρi isa Float64 &&
            default_info.min_floe_area isa Float64 &&
            default_info.min_floe_height isa Float64
        
        f32_info = FloeSettings(Float32)
        @test f32_info.ρi isa Float32 &&
            f32_info.min_floe_area isa Float32 &&
            f32_info.min_floe_height isa Float32
    end

    @testset "CouplingSettings" begin
        # Test default
        default_info = CouplingSettings()
        @test default_info.coupling_on &&
            default_info.Δt == 10 &&
            !default_info.two_way_coupling_on
        # Test custom correct
        custom_info = CouplingSettings(
            coupling_on = false,
            Δt = 20,
            two_way_coupling_on = false,
        )
        @test !custom_info.coupling_on &&
            custom_info.Δt == 20 &&
            !custom_info.two_way_coupling_on
        # Test negative timestep -> turns coupling off so turn calc_ocnτ off 
        warn_str1 = "Coupling can't occur on a multiple of negative timesteps. Turning coupling off."
        warn_str2 = "Can't calculate stresses on ocean from ice and atmosphere without coupling. Turning two_way_coupling_on off."
        neg_Δt_info = @test_logs (:warn, warn_str1) (:warn, warn_str2) CouplingSettings(
                coupling_on = true,
                Δt = -20,
                two_way_coupling_on = true,
            )
        @test !neg_Δt_info.coupling_on &&
            neg_Δt_info.Δt == -20 &&
            !neg_Δt_info.two_way_coupling_on
    end

    @testset "CollisionSettings" begin
        # Test default
        default_info = CollisionSettings()
        @test default_info.collisions_on &&
            default_info.floe_floe_max_overlap == 0.55 &&
            default_info.floe_domain_max_overlap == 0.75
        # Test custom false
        custom_info = CollisionSettings(
            collisions_on = false,
            floe_floe_max_overlap = -0.1,
            floe_domain_max_overlap = 2.0
        )
        @test !custom_info.collisions_on &&
            custom_info.floe_floe_max_overlap == -0.1 &&
            custom_info.floe_domain_max_overlap == 2
        # Test overlaps out of bounds
        f2large_str = "The maximum collisin overlap between floes can't be greater than 1. Setting to 1."
        f2small_str = "The maximum collisin overlap between floes can't be less than 0. Setting to 0."
        d2large_str = "The maximum collisin overlap between floes and the domain can't be greater than 1. Setting to 1."
        d2small_str = "The maximum collisin overlap between floes and the domain can't be less than 0. Setting to 0."
        out_of_bounds = @test_logs (:warn, f2large_str) (:warn, d2small_str) CollisionSettings(
                collisions_on = true,
                floe_floe_max_overlap = 1.2,
                floe_domain_max_overlap = -0.1,
            )
        @test out_of_bounds.collisions_on &&
            out_of_bounds.floe_floe_max_overlap == 1.0 &&
            out_of_bounds.floe_domain_max_overlap == 0.0
        
        out_of_bounds = @test_logs (:warn, f2small_str) (:warn, d2large_str) CollisionSettings(
                collisions_on = true,
                floe_floe_max_overlap = -0.1,
                floe_domain_max_overlap = 1.2,
            )
        @test out_of_bounds.collisions_on &&
            out_of_bounds.floe_floe_max_overlap == 0.0 &&
            out_of_bounds.floe_domain_max_overlap == 1.0

        @test typeof(CollisionSettings(Float32)) <: CollisionSettings{Float32}
        @test typeof(CollisionSettings(
            floe_floe_max_overlap = 1,
        )) <: CollisionSettings{Float64}

    end
    @testset "FractureSettings" begin
        # Test default
        default_info = FractureSettings()
        @test !default_info.fractures_on &&
            default_info.criteria isa NoFracture &&
            default_info.Δt == 0 &&
            !default_info.deform_on &&
            default_info.npieces == 3
        # Test custom correct
        criteria = HiblerYieldCurve(
            0.0,
            0.0,
            Subzero.make_polygon([[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]]),
        )
        custom_info = FractureSettings(
            fractures_on = true,
            criteria = criteria,
            Δt = 200,
            deform_on = true,
            npieces = 4
        )
        @test custom_info.fractures_on &&
            custom_info.criteria isa HiblerYieldCurve &&
            custom_info.Δt == 200 &&
            custom_info.deform_on &&
            custom_info.npieces == 4
        # Test warnings
        negΔt_str = "Fracturing can't occur with negative timesteps. Turning \
            fracturing off."
        nocriteria_str = "Fracturing can't occur on with NoFracture criteria. \
            Turning fracturing off."
        deform_str= "Deformation can't occur on without fracturing. Turning \
            deformation off."
        npiece_str = "Fracturing can't occur on with npieces < 2 as this won't \
            split floe. Turning fracturing off."
        # Test negative timestep provided -> can't fracture
        negΔt_info = @test_logs (:warn, negΔt_str) (:warn, deform_str) FractureSettings(
                fractures_on = true,
                criteria = criteria,
                Δt = -200,
                deform_on = true,
            )
        @test !negΔt_info.fractures_on &&
            negΔt_info.Δt == -200 &&
            !negΔt_info.deform_on
        # Test no fracture criteria provided -> can't fracture
        nocriteria_info = @test_logs (:warn, nocriteria_str) (:warn, deform_str) FractureSettings(
                fractures_on = true,
                criteria = NoFracture(),
                Δt = 200,
                deform_on = true,
            )
        @test !nocriteria_info.fractures_on && !nocriteria_info.deform_on
        # Test npieces == 1 -> won't fracture
        small_npiece_info = @test_logs (:warn, npiece_str) (:warn, deform_str) FractureSettings(
            fractures_on = true,
            criteria = criteria,
            Δt = 200,
            deform_on = true,
            npieces = 1,
        )
        @test !small_npiece_info.fractures_on && !small_npiece_info.deform_on
    end

    @testset "SimplificationSettings" begin
        # Test defaults
        default_info = SimplificationSettings()
        @test default_info.smooth_vertices_on &&
            default_info.max_vertices == 30 &&
            default_info.Δt_smooth == 20 &&
            default_info.tol == 100
        
        neg_Δt_str = "Floe smoothing can't occur on a multiple of negative \
            timesteps. Turning floe simplification off."
        neg_Δt_info = @test_logs (:warn, neg_Δt_str) SimplificationSettings(Δt_smooth = -20)
        @test !neg_Δt_info.smooth_vertices_on

        @test typeof(SimplificationSettings(Float32)) <: SimplificationSettings{Float32}
        @test typeof(SimplificationSettings(
            tol = 50,
        )) <: SimplificationSettings{Float64}
    end

    @testset "RidgeRaftSettings" begin
        default_info = RidgeRaftSettings()
        @test !default_info.ridge_raft_on &&
            default_info.Δt == 0 &&
            default_info.ridge_probability == 0.95 &&
            default_info.raft_probability == 0.95 &&
            default_info.min_overlap_frac == 0.01 &&
            default_info.min_ridge_height == 0.2 &&
            default_info.max_floe_ridge_height == 5.0 &&
            default_info.max_domain_ridge_height == 1.25 &&
            default_info.max_floe_raft_height == 0.25 &&
            default_info.max_domain_raft_height == 0.25 &&
            default_info.domain_gain_probability == 1.0
            
        neg_Δt_str = "Ridging and rafting can't occur on a multiple of negative \
            timesteps. Turning ridging and rafting off."
        neg_Δt_info = @test_logs (:warn, neg_Δt_str) RidgeRaftSettings(ridge_raft_on = true, Δt = -10)
        @test !neg_Δt_info.ridge_raft_on

        ridge_prob2_str = "Floes can't have a greater ridge probability than 1. \
            Setting ridge probability to 1."
        ridge_prob2_info = @test_logs (:warn, ridge_prob2_str) RidgeRaftSettings(ridge_raft_on = true, Δt = 100, ridge_probability = 2.0)
        @test ridge_prob2_info.ridge_probability == 1.0
    end
    @testset "WeldSettings" begin
        default_info = WeldSettings()
        @test !default_info.weld_on &&
            isempty(default_info.Δts) &&
            isempty(default_info.Nxs) &&
            isempty(default_info.Nys) &&
            default_info.min_weld_area == 1e6 &&
            default_info.max_weld_area == 2e9 &&
            default_info.welding_coeff == 150

        no_Δts_str = "Welding can't occur without any given timesteps or with \
            negative timesteps. Turning welding off."
        @test_logs (:warn, no_Δts_str) WeldSettings(weld_on = true)
        zero_split_str = "Can't split the grid into less than one row or column. \
            Turning welding off."
        @test_logs (:warn, zero_split_str) WeldSettings(
            weld_on = true, Δts = [100],
            Nxs = [1], Nys = [0],
        )
        length_wrong_str = "Length of timestep multiple list (Δts) must match \
        length of grid split lists Nxs and Nys. Turning welding off." 
        @test_logs (:warn, length_wrong_str) WeldSettings(
            weld_on = true, Δts = [100],
            Nxs = [1], Nys = [1, 2],
        )
        # Test sorting of timesteps
        sorted_settings = WeldSettings(
            weld_on = true,
            Δts = [100, 400, 700],
            Nxs = [1, 2, 3],
            Nys = [4, 5, 6],
        )
        @test sorted_settings.Δts == [700, 400, 100]
        @test sorted_settings.Nxs == [3, 2, 1]
        @test sorted_settings.Nys == [6, 5, 4]

        @test typeof(WeldSettings(Float32)) <: WeldSettings{Float32}
        @test typeof(WeldSettings(
            welding_coeff = 100,
        )) <: WeldSettings{Float64}
    end
end