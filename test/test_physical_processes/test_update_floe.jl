#=
Reference implementation of `timestep_floe_properties!`, copied from the CPU version on
`main` (including `calc_stress!`, `calc_strain!`, and `_move_floe!`) so that GPU/kernel
versions can be checked against it. Only supports `DecayAreaScaledCalculator`.
=#
function _reference_calc_stress!(floe, floe_settings, ::Type{FT}) where FT
    xi, yi = floe.centroid
    inters = floe.interactions
    stress = zeros(FT, 2, 2)
    if floe.num_inters > 0
        for i in 1:floe.num_inters
            stress[1, 1] += (inters[i, xpoint] - xi) * inters[i, xforce]
            stress[1, 2] += (inters[i, ypoint] - yi) * inters[i, xforce] +
                (inters[i, xpoint] - xi) * inters[i, yforce]
            stress[2, 2] += (inters[i, ypoint] - yi) * inters[i, yforce]
        end
        stress[1, 2] *= FT(0.5)
        stress[2, 1] = stress[1, 2]
        stress .*= 1/(floe.area * floe.height)
    end
    λ = floe_settings.stress_calculator.λ
    floe.stress_accum .= (1 - λ) * floe.stress_accum + λ * stress
    floe.stress_instant .= stress
    return
end

function _reference_calc_strain!(floe, ::Type{FT}) where FT
    fill!(floe.strain, zero(FT))
    trans_poly = Subzero._translate_poly(FT, floe.poly, -floe.centroid[1], -floe.centroid[2])
    local x1, y1
    for (i, p2) in enumerate(GI.getpoint(GI.getexterior(trans_poly)))
        x2, y2 = GO._tuple_point(p2, FT)
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        xdiff, ydiff = x2 - x1, y2 - y1
        rad1, rad2 = sqrt(x1^2 + y1^2), sqrt(x2^2 + y2^2)
        θ1, θ2 = atan(y1, x1), atan(y2, x2)
        u1 = floe.u - floe.ξ * rad1 * sin(θ1)
        u2 = floe.u - floe.ξ * rad2 * sin(θ2)
        v1 = floe.u + floe.ξ * rad1 * cos(θ1)  # sic: uses u, as on main
        v2 = floe.u + floe.ξ * rad2 * cos(θ2)
        udiff, vdiff = u2 - u1, v2 - v1
        floe.strain[1, 1] += udiff * ydiff
        floe.strain[1, 2] += udiff * xdiff + vdiff * ydiff
        floe.strain[2, 2] += vdiff * xdiff
        x1, y1 = x2, y2
    end
    floe.strain[1, 2] *= FT(0.5)
    floe.strain[2, 1] = floe.strain[1, 2]
    floe.strain ./= 2floe.area
    return
end

function _reference_move_floe!(floe, Δx, Δy, Δα, ::Type{FT}) where FT
    CT = Subzero.CoordinateTransformations
    cx, cy = floe.centroid
    floe.centroid[1] += Δx
    floe.centroid[2] += Δy
    rot = CT.LinearMap(Subzero.Rotations.Angle2d(Δα))
    cent_rot = CT.recenter(rot, (cx, cy))
    trans = CT.Translation(Δx, Δy)
    floe.poly = GO.tuples(GO.transform(trans ∘ cent_rot, floe.poly), FT)
    return
end

function _reference_timestep_floe_properties!(
    floes::StructArray{<:Floe{FT}},
    tstep,
    Δt,
    floe_settings,
) where FT
    for i in eachindex(floes)
        cforce = floes.collision_force[i]
        ctrq = floes.collision_trq[i]
        _reference_calc_stress!(Subzero.get_floe(floes, i), floe_settings, FT)

        if floes.height[i] > floe_settings.max_floe_height
            @info "Reducing height to $(floe_settings.max_floe_height) m" tstep = tstep
            floes.height[i] = floe_settings.max_floe_height
        end

        while maximum(abs.(cforce)) > floes.mass[i]/(5Δt)
            @info "Decreasing collision forces by a factor of 10" tstep = tstep
            cforce = cforce ./ 10
            ctrq = ctrq ./ 10
        end

        h = floes.height[i]
        Δh = floes.hflx_factor[i] / h
        hfrac = (h + Δh) / h
        floes.mass[i] *= hfrac
        floes.moment[i] *= hfrac
        floes.height[i] -= Δh
        h = floes.height[i]

        Δx = 1.5Δt*floes.u[i] - 0.5Δt*floes.p_dxdt[i]
        Δy = 1.5Δt*floes.v[i] - 0.5Δt*floes.p_dydt[i]
        Δα = 1.5Δt*floes.ξ[i] - 0.5Δt*floes.p_dαdt[i]
        floes.α[i] += Δα

        _reference_move_floe!(Subzero.get_floe(floes, i), Δx, Δy, Δα, FT)
        floes.p_dxdt[i] = floes.u[i]
        floes.p_dydt[i] = floes.v[i]
        floes.p_dαdt[i] = floes.ξ[i]

        dudt = (floes.fxOA[i] + cforce[1])/floes.mass[i]
        dvdt = (floes.fyOA[i] + cforce[2])/floes.mass[i]
        frac = if abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) > (h/2)
            frac1 = (sign(dudt)*h/2Δt)/dudt
            frac2 = (sign(dvdt)*h/2Δt)/dvdt
            min(frac1, frac2)
        elseif abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) < (h/2)
            (sign(dudt)*h/2Δt)/dudt
        elseif abs(Δt*dudt) < (h/2) && abs(Δt*dvdt) > (h/2)
            (sign(dvdt)*h/2Δt)/dvdt
        else
            1
        end
        if frac != 1
            @info "Adjusting u and v velocities to prevent too high" tstep = tstep
            dudt = frac*dudt
            dvdt = frac*dvdt
        end
        floes.u[i] += 1.5Δt*dudt-0.5Δt*floes.p_dudt[i]
        floes.v[i] += 1.5Δt*dvdt-0.5Δt*floes.p_dvdt[i]
        floes.p_dudt[i] = dudt
        floes.p_dvdt[i] = dvdt

        dξdt = (floes.trqOA[i] + ctrq)/floes.moment[i]
        dξdt = frac*dξdt
        ξ = floes.ξ[i] + 1.5Δt*dξdt-0.5Δt*floes.p_dξdt[i]
        if abs(ξ) > floe_settings.maximum_ξ
            @info "Shrinking ξ" tstep = tstep
            ξ = sign(ξ) * floe_settings.maximum_ξ
        end
        floes.ξ[i] = ξ
        floes.p_dξdt[i] = dξdt

        _reference_calc_strain!(Subzero.get_floe(floes, i), FT)
    end
    return
end

@testset "Update floe" begin
    backends = Any[Subzero.KernelAbstractions.CPU()]
    CUDA.functional() && push!(backends, CUDA.CUDABackend())
    @testset "timestep_floe_properties on $backend" for backend in backends
        FT = Float64
        Δt, tstep = 10, 1
        floe_settings = FloeSettings()
        expected = _make_timestep_test_floes(FT, floe_settings)
        floes = deepcopy(expected)
        with_logger(NullLogger()) do
            _reference_timestep_floe_properties!(expected, tstep, Δt, floe_settings)
        end
        @test_logs(
            (:info, r"Reducing height"),
            (:info, "Decreasing collision forces by a factor of 10"),
            (:info, "Adjusting u and v velocities to prevent too high"),
            (:info, "Shrinking ξ"),
            match_mode = :any,
            Subzero.timestep_floe_properties!(floes, tstep, Δt, floe_settings; backend),
        )
        rtol = 1e-10
        @testset "floe $i" for i in eachindex(floes)
            points = collect(GI.getpoint(floes.poly[i]))
            expected_points = collect(GI.getpoint(expected.poly[i]))
            @test length(points) == length(expected_points)
            if length(points) == length(expected_points)
                @test isapprox(
                    reinterpret(FT, points),
                    reinterpret(FT, expected_points);
                    rtol,
                )
            end
            @testset "$field" for field in (
                :centroid, :height, :mass, :moment, :α, :u, :v, :ξ,
                :p_dxdt, :p_dydt, :p_dαdt, :p_dudt, :p_dvdt, :p_dξdt,
                :stress_accum, :stress_instant, :strain,
            )
                actual_value = getproperty(floes, field)[i]
                expected_value = getproperty(expected, field)[i]
                @test isapprox(actual_value, expected_value; rtol)
            end
        end
    end
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
                floe_dict["height"][i];
                u = floe_dict["u"][i],
                v = floe_dict["v"][i],
                ξ = floe_dict["ξ"][i],
                floe_settings = floe_settings,
            )
            f.interactions = floe_dict["interactions"][i]
            f.num_inters = size(f.interactions, 1)
            f.stress_instant = floe_dict["last_stress"][i]
            stress = Subzero.calc_stress!(f, floe_settings)
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
            @test f.poly == Subzero.make_polygon(floe_dict["coords"][i])
        end
    end
    @testset "calc_stress! on FixedWidthFloes" begin
        for FT in (Float64, Float32)
            floe_settings = FloeSettings(FT)
            floes = _make_timestep_test_floes(FT, floe_settings)
            fwf = Subzero.FixedWidthFloes(floes)
            for i in eachindex(floes)
                Subzero.calc_stress!(Subzero.get_floe(floes, i), floe_settings)
                Subzero.calc_stress!(fwf, i, floe_settings.stress_calculator)
                @test isapprox(fwf.stress_instant[i, :, :], floes.stress_instant[i]; rtol = 10eps(FT))
                @test isapprox(fwf.stress_accum[i, :, :], floes.stress_accum[i]; rtol = 10eps(FT))
            end
        end
    end
    @testset "calc_strain! on FixedWidthFloes" begin
        for FT in (Float64, Float32)
            floes = _make_timestep_test_floes(FT, FloeSettings(FT))
            fwf = Subzero.FixedWidthFloes(floes)
            for i in eachindex(floes)
                Subzero.calc_strain!(Subzero.get_floe(floes, i))
                Subzero.calc_strain!(fwf, i)
                @test isapprox(fwf.strain[i, :, :], floes.strain[i]; rtol = 10eps(FT))
            end
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
        f1 = Floe(coords1, 0.5)  # this is a square
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
        @test GO.equals(f1.poly, tri_poly)
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
            0.5;
            u = 0.1,
            v = 0.25,
            ξ = -0.5,
        )
        sqr_floe.p_dxdt = 0.11
        sqr_floe.p_dydt = 0.22
        sqr_floe.p_dαdt = -0.45

        tri_floe = Floe(
            triangle_coords,
            0.5;
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
        translate_coords!(triangle_coords, 10.0, 0.0)
        sqr_floe = Floe(
            square_coords,
            0.5;
            u = 0.1,
            v = 0.25,
            ξ = -0.5,
        )
        sqr_floe.p_dxdt = 0.11
        sqr_floe.p_dydt = 0.22
        sqr_floe.p_dαdt = -0.45

        tri_floe = Floe(
            triangle_coords,
            0.5;
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
            0.5
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
            Floe(left_coords, 0.5),
            Floe(right_and_mid_coords, 0.5)
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
            Floe(left_coords, 0.5),
            Floe(right_coords, 0.5),
            Floe(mid_coords, 0.5)
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