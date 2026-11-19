#=
The collision forces between two floes on the device (the `FixedWidthFloes` methods in
`src/physical_processes/collisions.jl`) must be the same as `floe_floe_interaction!`, bit for
bit. These tests compare them pair by pair on the CPU, both with the CPU functions that
`floe_floe_interaction!` calls (which also gives the forces that are zero) and with
`floe_floe_interaction!` itself, and then run them in a kernel on every backend.
=#
@testset "Device collision forces" begin
    Δt = 10
    #= In the floes' type, as in a simulation (`CollisionSettings{FT}`): the CPU
    `floe_floe_interaction!` takes its float type from `max_overlap` =#
    _max_overlap(FT) = FT(CollisionSettings().floe_floe_max_overlap)

    #= Pairs of floes that overlap by different amounts: real floe shapes next to each other
    (by how much along x is random, so a few overlap enough to fuse), and some synthetic
    polygons, which mostly overlap a lot. Floes 2k - 1 and 2k are a pair. The floes move and
    turn, so that there is friction. =#
    function _collision_test_floes(::Type{FT}; rng = Xoshiro(21)) where FT
        fs = FloeSettings(FT; min_floe_area = 0)
        shapes = _floe_shape_rings()
        # every sixth shape with up to 300 vertices, to keep the buffers small
        shapes = filter(s -> length(s) <= 300, shapes)[1:6:end]
        polys = Vector{Vector{Tuple{Float64, Float64}}}()
        # for each shape, a pair that overlaps a little and one that overlaps more
        for p in shapes, depth in (0.6, 0.3)
            q = rand(rng, shapes)
            cp, cq = (mean(first, p), mean(last, p)), (mean(first, q), mean(last, q))
            half_widths = (maximum(first, p) - minimum(first, p) +
                maximum(first, q) - minimum(first, q)) / 2
            Δx = cp[1] - cq[1] + rand(rng, (-1, 1)) * half_widths * (depth + 0.4 * rand(rng))
            Δy = cp[2] - cq[2] + 0.1 * half_widths * randn(rng)
            push!(polys, p, [(x + Δx, y + Δy) for (x, y) in q])
        end
        append!(polys, _geom_test_polys(Float64; n = 20, rng))
        floes = StructArray([Floe{FT}([[[FT(x), FT(y)] for (x, y) in [p; [p[1]]]]],
            FT(0.25 + rand(rng)); floe_settings = fs, rng) for p in polys])
        floes.u .= FT.(0.5 .* randn(rng, length(floes)))
        floes.v .= FT.(0.5 .* randn(rng, length(floes)))
        floes.ξ .= FT.(1e-4 .* randn(rng, length(floes)))
        return floes
    end

    #= The forces between floes i and j from the CPU functions that `floe_floe_interaction!`
    calls, in the same way: whether they fuse, and the force, point and overlap of each
    region (rows), including forces that are zero =#
    function _cpu_pair_forces(floes, i, j, consts, max_overlap)
        FT = eltype(floes.area)
        ifloe, jfloe = Subzero.get_floe(floes, i), Subzero.get_floe(floes, j)
        regions = Subzero.intersect_polys(ifloe.poly, jfloe.poly, FT)
        areas = [FT(Subzero.area_poly(r, FT)) for r in regions]
        total_area = sum(areas; init = FT(0))
        empty_rows = (zeros(FT, 0, 2), zeros(FT, 0, 2), zeros(FT, 0))
        total_area > 0 || return false, empty_rows...
        max(total_area / ifloe.area, total_area / jfloe.area) > max_overlap &&
            return true, empty_rows...
        ih, jh = ifloe.height, jfloe.height
        ir, jr = FT(sqrt(ifloe.area)), FT(sqrt(jfloe.area))
        force_factor = ir > 1e5 || jr > 1e5 ? consts.E * min(ih, jh) / min(ir, jr) :
            consts.E * (ih * jh) / (ih * jr + jh * ir)
        normal, points, overlaps, Δl = Subzero.calc_elastic_forces(ifloe.poly, jfloe.poly,
            regions, areas, force_factor)
        size(points, 1) == 0 && return false, empty_rows...
        friction = Subzero.calc_friction_forces(ifloe, jfloe, points, normal, Δl, consts, Δt)
        return false, normal .+ friction, points, overlaps
    end

    # The forces of row k of `out` as the rows of _cpu_pair_forces
    _device_rows(out, k, n) = (out.force[k, 1:n, :], out.point[k, 1:n, :], out.overlap[k, 1:n])

    # Pairs both ways round, since floe i gets the forces
    _test_pairs(n) = [p for k in 1:2:(n - 1) for p in ((k, k + 1), (k + 1, k))]

    _scratch(floes, fwf, npairs) = (Subzero._CollisionScratch(eltype(floes.area), npairs,
        size(fwf.poly, 2), 4 * size(fwf.poly, 2) + 64, 8),
        Subzero._PairForces(eltype(floes.area), npairs, 8))

    @testset "Pairs of floes with $FT" for FT in (Float64, Float32)
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = _test_pairs(length(floes))
        scratch, out = _scratch(floes, fwf, length(pairs))
        same_rows, same_interactions = true, true
        n_fuse, n_forces, n_regions_many_points = 0, 0, 0
        for (k, (i, j)) in enumerate(pairs)
            n, fuse = Subzero.floe_floe_interaction!(Subzero._PairForces(out, k),
                Subzero._CollisionScratch(scratch, k), fwf, i, j, consts, Δt, max_overlap)
            @test n >= 0
            # with the CPU functions
            cpu_fuse, cpu_rows... = _cpu_pair_forces(floes, i, j, consts, max_overlap)
            same_rows &= fuse == cpu_fuse && isequal(_device_rows(out, k, n), cpu_rows)
            # with floe_floe_interaction! itself, which leaves out forces that are zero
            cpu_floes = deepcopy(floes[[i, j]])
            Subzero.floe_floe_interaction!(Subzero.get_floe(cpu_floes, 1), 1, Subzero.get_floe(cpu_floes, 2),
                2, consts, Δt, max_overlap)
            ninters = cpu_floes.num_inters[1]
            inters = cpu_floes.interactions[1][1:ninters, :]
            kept = [r for r in 1:n if out.force[k, r, 1] != 0 || out.force[k, r, 2] != 0]
            same_interactions &= (cpu_floes.status[1].tag == Subzero.fuse) == fuse &&
                isequal(inters[:, Subzero.xforce:Subzero.yforce], out.force[k, kept, :]) &&
                isequal(inters[:, Subzero.xpoint:Subzero.ypoint], out.point[k, kept, :]) &&
                isequal(inters[:, Subzero.overlap], out.overlap[k, kept])
            n_fuse += fuse
            n_forces += n
            # regions where the force direction comes from _many_intersect_normal_force!
            if n > 0
                ipoints = Subzero.get_intersection_points(floes.poly[i], floes.poly[j])
                for r in Subzero.intersect_polys(floes.poly[i], floes.poly[j], FT)
                    m = length(Subzero.which_vertices_match_points(ipoints, r))
                    n_regions_many_points += m != 0 && m != 2
                end
            end
        end
        @test same_rows
        @test same_interactions
        @test n_fuse > 3
        @test n_forces > 50
        @test n_regions_many_points > 10
    end

    @testset "Too little room with $FT" for FT in (Float64, Float32)
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        # the pair with the most forces
        pairs = _test_pairs(length(floes))
        scratch, out = _scratch(floes, fwf, length(pairs))
        counts = map(enumerate(pairs)) do (k, (i, j))
            first(Subzero.floe_floe_interaction!(Subzero._PairForces(out, k),
                Subzero._CollisionScratch(scratch, k), fwf, i, j, consts, Δt, max_overlap))
        end
        i, j = pairs[argmax(counts)]
        n = maximum(counts)
        @test n > 1
        clip_run(scratch, out) = first(Subzero.floe_floe_interaction!(Subzero._PairForces(out, 1),
            Subzero._CollisionScratch(scratch, 1), fwf, i, j, consts, Δt, max_overlap))
        big = Subzero._CollisionScratch(FT, 1, size(fwf.poly, 2), 4 * size(fwf.poly, 2) + 64, 8)
        @test clip_run(big, Subzero._PairForces(FT, 1, n)) == n
        @test clip_run(big, Subzero._PairForces(FT, 1, n - 1)) == Subzero._CLIP_TOO_SMALL
        small_lists = Subzero._CollisionScratch(FT, 1, size(fwf.poly, 2), 8, 8)
        @test clip_run(small_lists, Subzero._PairForces(FT, 1, n)) == Subzero._CLIP_TOO_SMALL
        # room for the floe's polygon but not for its moved copy
        no_room_to_move = Subzero._CollisionScratch(FT, 1, 4, 4 * size(fwf.poly, 2) + 64, 8)
        @test clip_run(no_room_to_move, Subzero._PairForces(FT, 1, n)) == Subzero._CLIP_TOO_SMALL
    end

    #= Kernel test: one thread per pair. `launch_per_floe!` only needs `n_floes`, so it gets a
    NamedTuple with the number of pairs instead of floes. =#
    function _pair_forces_kernel!(_, k, floes, pairs, scratch, out, results, consts, Δt,
            max_overlap)
        n, fuse = Subzero.floe_floe_interaction!(Subzero._PairForces(out, k),
            Subzero._CollisionScratch(scratch, k), floes, pairs[k, 1], pairs[k, 2], consts,
            Δt, max_overlap)
        results[k, 1] = n
        results[k, 2] = fuse
        return
    end
    @testset "Pairs of floes on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = Int32.(reduce(vcat, [[i j] for (i, j) in _test_pairs(length(floes))]))
        npairs = size(pairs, 1)
        outputs() = (_scratch(floes, fwf, npairs)..., zeros(Int32, npairs, 2))
        ref = outputs()
        for k in 1:npairs
            _pair_forces_kernel!(nothing, k, fwf, pairs, ref..., consts, Δt, max_overlap)
        end
        @test all(>=(0), ref[3][:, 1])
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_pair_forces_kernel!, backend, (; n_floes = npairs),
            Subzero.adapt(backend, fwf), Subzero.adapt(backend, pairs), dev..., consts, Δt,
            max_overlap)
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test got[3] == ref[3]
        @test all(1:npairs) do k
            n = ref[3][k, 1]
            isequal(_device_rows(got[2], k, n), _device_rows(ref[2], k, n))
        end
        @test sum(ref[3][:, 1]) > 50
    end
end
