@testset "Update floe" begin
    @testset "Stress/Strain" begin
        # Test Stress and Strain Calculations
        floe_dict = load(
            "inputs/stress_strain.jld2"  # uses the first 2 element
        )
        floe_settings = FloeSettings()
        stresses = [[-10.065, 36.171, 36.171, -117.458],
            [7.905, 21.913, 21.913, -422.242]]
        stress_histories = [[-4971.252, 17483.052, 17483.052, -57097.458],
            [4028.520, 9502.886, 9502.886, -205199.791]]
        strains = [[-0.0372, 0, 0, .9310], [7.419, 0, 0, -6.987]]
        strain_multiplier = [1e6, 1e6]
        for i in 1:2
            f = Floe(
                floe_dict["coords"][i],
                floe_dict["height"][i],
                0.0,
                u = floe_dict["u"][i],
                v = floe_dict["v"][i],
                ξ = floe_dict["ξ"][i],
                floe_settings = floe_settings,
            )
            f.interactions = floe_dict["interactions"][i]
            f.num_inters = size(f.interactions, 1)
            f.stress_instant = floe_dict["last_stress"][i]
            stress = Subzero.calc_stress!(f, floe_settings, 10)
            @test_broken all(isapprox.(vec(f.stress_accum), stresses[i], atol = 1e-3))
            @test all(isapprox.(
                vec(f.stress_instant),
                stress_histories[i],
                atol = 1e-3
            ))
            Subzero.calc_strain!(f)
            @test all(isapprox.(
                vec(f.strain) .* strain_multiplier[i],
                strains[i],
                atol = 1e-3
            ))
            @test f.coords == floe_dict["coords"][i]
        end
    end
    @testset "Replace floe" begin
        # Test replace floe
        coords1 = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
            [0.0, 0.0],
        ]]
        f1 = Floe(coords1, 0.5, 0.0)  # this is a square
        mass1 = f1.mass
        triangle_coords = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [0.0, 0.0],
        ]]
        tri_poly = Subzero.make_polygon(triangle_coords)  # this is a triangle
        Subzero.replace_floe!(
            f1,
            tri_poly,
            f1.mass,
            FloeSettings(),
            Xoshiro(1)
        )
        @test all(f1.centroid .== GO.centroid(tri_poly))
        @test f1.coords == Subzero.find_poly_coords(tri_poly)
        @test f1.area == GO.area(tri_poly)
        @test f1.mass == mass1
        @test f1.height * f1.area * 920.0 == f1.mass
        @test f1.α == 0
        @test f1.status.tag == Subzero.active
        @test f1.rmax == 10*√5 / 3
    end
    @testset "Conserve momentum" begin
        square_coords = [[
            [0.0, 0.0],
            [0.0, 20.0],
            [20.0, 20.0],
            [20.0, 0.0],
            [0.0, 0.0],
        ]]
        triangle_coords = [[
            [0.0, 0.0],
            [10.0, 20.0],
            [20.0, 0.0],
            [0.0, 0.0],
        ]]
        sqr_floe = Floe(
            square_coords,
            0.5,
            0.0,
            u = 0.1,
            v = 0.25,
            ξ = -0.5,
        )
        sqr_floe.p_dxdt = 0.11
        sqr_floe.p_dydt = 0.22
        sqr_floe.p_dαdt = -0.45

        tri_floe = Floe(
            triangle_coords,
            0.5,
            0.0,
            u = 0.1,
            v = 0.25,
            ξ = -0.5,
        )
        tri_floe.p_dxdt = 0.11
        tri_floe.p_dydt = 0.22
        tri_floe.p_dαdt = -0.45
        # Test one floe changing shape
        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            [sqr_floe.u],
            [sqr_floe.v],
            [sqr_floe.mass],
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            [sqr_floe.u],
            [sqr_floe.v],
            [sqr_floe.mass],
            [sqr_floe.ξ],
            [sqr_floe.moment],
            [sqr_floe.centroid[1]],
            [sqr_floe.centroid[2]],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            [sqr_floe.p_dxdt],
            [sqr_floe.p_dydt],
            [sqr_floe.mass],
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            [sqr_floe.p_dxdt],
            [sqr_floe.p_dydt],
            [sqr_floe.mass],
            [sqr_floe.p_dαdt],
            [sqr_floe.moment],
            [sqr_floe.centroid[1]],
            [sqr_floe.centroid[2]],
        )
        Subzero.conserve_momentum_change_floe_shape!(
            sqr_floe.mass,
            sqr_floe.moment,
            sqr_floe.centroid[1],
            sqr_floe.centroid[2],
            10,
            tri_floe,
        )
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            [tri_floe.u],
            [tri_floe.v],
            [tri_floe.mass],
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            [tri_floe.u],
            [tri_floe.v],
            [tri_floe.mass],
            [tri_floe.ξ],
            [tri_floe.moment],
            [tri_floe.centroid[1]],
            [tri_floe.centroid[2]],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            [tri_floe.p_dxdt],
            [tri_floe.p_dydt],
            [tri_floe.mass],
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            [tri_floe.p_dxdt],
            [tri_floe.p_dydt],
            [tri_floe.mass],
            [tri_floe.p_dαdt],
            [tri_floe.moment],
            [tri_floe.centroid[1]],
            [tri_floe.centroid[2]],
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-8)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-8)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-8)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-8)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-8,
        )
        @test isapprox(
            p_spin_momentum_init + p_angular_momentum_init,
            p_spin_momentum_after + p_angular_momentum_after,
            atol = 1e-8,
        )
        # Test two floes combining
        Subzero.translate!(triangle_coords, 10.0, 0.0)
        sqr_floe = Floe(
            square_coords,
            0.5,
            0.0,
            u = 0.1,
            v = 0.25,
            ξ = -0.5,
        )
        sqr_floe.p_dxdt = 0.11
        sqr_floe.p_dydt = 0.22
        sqr_floe.p_dαdt = -0.45

        tri_floe = Floe(
            triangle_coords,
            0.5,
            0.0,
            u = 0.3,
            v = 0.05,
            ξ = 0.2,
        )
        tri_floe.p_dxdt = 0.2
        tri_floe.p_dydt = 0.04
        tri_floe.p_dαdt = 0.19

        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            [sqr_floe.u, tri_floe.u],
            [sqr_floe.v, tri_floe.v],
            [sqr_floe.mass, tri_floe.mass],
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            [sqr_floe.u, tri_floe.u],
            [sqr_floe.v, tri_floe.v],
            [sqr_floe.mass, tri_floe.mass],
            [sqr_floe.ξ, tri_floe.ξ],
            [sqr_floe.moment, tri_floe.moment],
            [sqr_floe.centroid[1], tri_floe.centroid[1]],
            [sqr_floe.centroid[2], tri_floe.centroid[2]],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            [sqr_floe.p_dxdt, tri_floe.p_dxdt],
            [sqr_floe.p_dydt, tri_floe.p_dydt],
            [sqr_floe.mass, tri_floe.mass],
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            [sqr_floe.p_dxdt, tri_floe.p_dxdt],
            [sqr_floe.p_dydt, tri_floe.p_dydt],
            [sqr_floe.mass, tri_floe.mass],
            [sqr_floe.p_dαdt, tri_floe.p_dαdt],
            [sqr_floe.moment, tri_floe.moment],
            [sqr_floe.centroid[1], tri_floe.centroid[1]],
            [sqr_floe.centroid[2], tri_floe.centroid[2]],
        )
        mass1 = sqr_floe.mass
        moment1 = sqr_floe.moment
        x1, y1 = sqr_floe.centroid
        tri_rect_poly = Subzero.union_polys(Subzero.make_polygon(square_coords), Subzero.make_polygon(triangle_coords))[1]
        Subzero.replace_floe!(
            sqr_floe,
            tri_rect_poly,
            sqr_floe.mass + tri_floe.mass,
            FloeSettings(),
            Xoshiro(1)
        )
        Subzero.conserve_momentum_change_floe_shape!(
            mass1,
            moment1,
            x1,
            y1,
            10,
            sqr_floe,
            tri_floe,
        )
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            [sqr_floe.u],
            [sqr_floe.v],
            [sqr_floe.mass],
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            [sqr_floe.u],
            [sqr_floe.v],
            [sqr_floe.mass],
            [sqr_floe.ξ],
            [sqr_floe.moment],
            [sqr_floe.centroid[1]],
            [sqr_floe.centroid[2]],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            [sqr_floe.p_dxdt],
            [sqr_floe.p_dydt],
            [sqr_floe.mass],
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            [sqr_floe.p_dxdt],
            [sqr_floe.p_dydt],
            [sqr_floe.mass],
            [sqr_floe.p_dαdt],
            [sqr_floe.moment],
            [sqr_floe.centroid[1]],
            [sqr_floe.centroid[2]],
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-8)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-8)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-8)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-8)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-8,
        )
        @test isapprox(
            p_spin_momentum_init + p_angular_momentum_init,
            p_spin_momentum_after + p_angular_momentum_after,
            atol = 1e-8,
        )

        # Fracturing a floe momentum conservation
        initial_coords = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [20.0, 10.0],
            [20.0, 0.0],
            [0.0, 0.0],
        ]]
        left_coords = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [5.0, 10.0],
            [5.0,  0.0],
            [0.0, 0.0],
        ]]
        mid_coords = [[
            [5.0, 0.0],
            [5.0, 10.0],
            [15.0, 10.0],
            [15.0, 0.0],
            [5.0, 0.0],
        ]]
        right_coords = [[
            [15.0, 0.0],
            [15.0, 10.0],
            [20.0, 10.0],
            [20.0, 0.0],
            [15.0, 0.0],
        ]]
        right_and_mid_coords = [[
            [5.0, 0.0],
            [5.0, 10.0],
            [20.0, 10.0],
            [20.0, 0.0],
            [5.0, 0.0],
        ]]
        #  One floe splitting into two floes
        initial_floe = Floe(
            initial_coords,
            0.5,
            0.0
        )
        initial_floe.u = 0.1
        initial_floe.v = -0.2
        initial_floe.ξ = -0.08
        initial_floe.p_dαdt = -0.03
        initial_floe.p_dxdt = 0.09

        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            [initial_floe.u],
            [initial_floe.v,],
            [initial_floe.mass],
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            [initial_floe.u],
            [initial_floe.v,],
            [initial_floe.mass],
            [initial_floe.ξ],
            [initial_floe.moment],
            [initial_floe.centroid[1]],
            [initial_floe.centroid[2]],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            [initial_floe.p_dxdt],
            [initial_floe.p_dydt],
            [initial_floe.mass],
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            [initial_floe.p_dxdt],
            [initial_floe.p_dydt],
            [initial_floe.mass],
            [initial_floe.p_dαdt],
            [initial_floe.moment],
            [initial_floe.centroid[1]],
            [initial_floe.centroid[2]],
        )
        new_floes = StructArray([
            Floe(
                left_coords,
                0.5,
                0.0
            ),
            Floe(
                right_and_mid_coords,
                0.5,
                0.0,
            )
        ])
        Subzero.conserve_momentum_fracture_floe!(
            initial_floe,
            new_floes,
            10,
        )
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            new_floes.u,
            new_floes.v,
            new_floes.mass,
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            new_floes.u,
            new_floes.v,
            new_floes.mass,
            new_floes.ξ,
            new_floes.moment,
            [c[1] for c in new_floes.centroid],
            [c[2] for c in new_floes.centroid],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            new_floes.p_dxdt,
            new_floes.p_dydt,
            new_floes.mass,
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            new_floes.p_dxdt,
            new_floes.p_dydt,
            new_floes.mass,
            new_floes.p_dαdt,
            new_floes.moment,
            [c[1] for c in new_floes.centroid],
            [c[2] for c in new_floes.centroid],
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-8)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-8)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-8)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-8)

        # One floe splitting into three floes
        new_floes = StructArray([
            Floe(
                left_coords,
                0.5,
                0.0
            ),
            Floe(
                right_coords,
                0.5,
                0.0
            ),
            Floe(
                mid_coords,
                0.5,
                0.0,
            )
        ])
        Subzero.conserve_momentum_fracture_floe!(
            initial_floe,
            new_floes,
            10,
        )
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            new_floes.u,
            new_floes.v,
            new_floes.mass,
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            new_floes.p_dxdt,
            new_floes.p_dydt,
            new_floes.mass,
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-8)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-8)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-8)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-8)
    end
end