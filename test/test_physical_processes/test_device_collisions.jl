#=
The collision forces between two floes on the device (the `floe_floe_interaction!` method for
a `FixedWidthFloes`, and the functions it calls) against the CPU functions that
`floe_floe_interaction!` calls for a StructArray of floes: bit for bit with Float64 floes,
and approximately with Float32 floes, for which the CPU functions compute the intersection
points, and some steps of the clipping, in Float64.
=#
@testset "Device collision forces" begin
    Δt = 10
    #= In the floes' type, as in a simulation (`CollisionSettings{FT}`): the CPU
    `floe_floe_interaction!` takes its float type from `max_overlap` =#
    _max_overlap(FT) = FT(CollisionSettings().floe_floe_max_overlap)

    #= Pairs of floes that overlap by different amounts: real floe shapes next to each other
    (by how much along x is random, so a few overlap enough to fuse), and some synthetic
    polygons, which mostly overlap a lot. Floes 2k - 1 and 2k are a pair. The floes move and
    turn, so that there is friction. The floes are the same for both float types, apart from
    rounding: making a floe uses random numbers for its sub-floe points, how many depending on
    the float type, so it gets its own generator. =#
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
        n = length(polys)
        heights, u, v, ξ = 0.25 .+ rand(rng, n), 0.5 .* randn(rng, n), 0.5 .* randn(rng, n),
            1e-4 .* randn(rng, n)
        floes = StructArray([Floe{FT}([[[FT(x), FT(y)] for (x, y) in [p; [p[1]]]]],
            FT(heights[k]); floe_settings = fs, rng = Xoshiro(k)) for (k, p) in enumerate(polys)])
        floes.u .= FT.(u)
        floes.v .= FT.(v)
        floes.ξ .= FT.(ξ)
        return floes
    end

    #= The forces on floe i from floe j from the CPU functions that `floe_floe_interaction!`
    calls, in the same way: whether they fuse, and the force, point and overlap of each region
    (rows), including forces that are zero =#
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

    # Room for the real floe shapes, which can intersect at many points
    _buffers(fwf, npairs; max_ipoints = 64, max_regions = 8) = (
        Subzero.CollisionWorkspace(eltype(fwf.area), npairs, size(fwf.poly, 2); max_ipoints,
            max_regions),
        Subzero.PairForces(eltype(fwf.area), npairs, max_regions))

    _pair_run(out, ws, fwf, k, i, j, consts, max_overlap) = Subzero.floe_floe_interaction!(
        Subzero.PairForces(out, k), Subzero.CollisionWorkspace(ws, k), fwf, i, j, consts, Δt,
        max_overlap)

    #= The overlap regions of floes i and j from the CPU functions, each with the number of
    intersection points m that `calc_normal_force` finds in it (0 if the floes' edges
    intersect at fewer than two points, when there are no forces). The force direction comes
    from `_many_intersect_normal_force` where m isn't 0 or 2. =#
    function _regions_and_points(floes, i, j)
        FT = eltype(floes.area)
        ipoints = Subzero.get_intersection_points(floes.poly[i], floes.poly[j])
        regions = Subzero.intersect_polys(floes.poly[i], floes.poly[j], FT)
        length(ipoints) >= 2 || return [(r, 0) for r in regions]
        return [(r, length(Subzero.which_vertices_match_points(ipoints, r))) for r in regions]
    end
    _many_points(m) = m != 0 && m != 2

    @testset "Pairs of floes with Float64" begin
        FT = Float64
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = _test_pairs(length(floes))
        ws, out = _buffers(fwf, length(pairs))
        # pairs whose forces differ from the CPU functions' and floe_floe_interaction!'s
        different_rows, different_interactions = Tuple{Int, Int}[], Tuple{Int, Int}[]
        n_fuse, n_forces, n_many = 0, 0, 0
        for (k, (i, j)) in enumerate(pairs)
            n, fuse = _pair_run(out, ws, fwf, k, i, j, consts, max_overlap)
            @test n >= 0
            # with the CPU functions
            cpu_fuse, cpu_rows... = _cpu_pair_forces(floes, i, j, consts, max_overlap)
            fuse == cpu_fuse && isequal(_device_rows(out, k, n), cpu_rows) ||
                push!(different_rows, (i, j))
            # with floe_floe_interaction! itself, which leaves out forces that are zero
            cpu_floes = deepcopy(floes[[i, j]])
            Subzero.floe_floe_interaction!(Subzero.get_floe(cpu_floes, 1), 1,
                Subzero.get_floe(cpu_floes, 2), 2, consts, Δt, max_overlap)
            ninters = cpu_floes.num_inters[1]
            inters = cpu_floes.interactions[1][1:ninters, :]
            kept = [r for r in 1:n if out.force[k, r, 1] != 0 || out.force[k, r, 2] != 0]
            (cpu_floes.status[1].tag == Subzero.fuse) == fuse &&
                isequal(inters[:, Subzero.xforce:Subzero.yforce], out.force[k, kept, :]) &&
                isequal(inters[:, Subzero.xpoint:Subzero.ypoint], out.point[k, kept, :]) &&
                isequal(inters[:, Subzero.overlap], out.overlap[k, kept]) ||
                push!(different_interactions, (i, j))
            n_fuse += fuse
            n_forces += n
            n > 0 && (n_many += count(((_, m),) -> _many_points(m), _regions_and_points(floes, i, j)))
        end
        @test isempty(different_rows)
        @test isempty(different_interactions)
        # The test covers floes that fuse, many forces, and regions with many intersections
        @test n_fuse > 3
        @test n_forces > 50
        @test n_many > 10
    end

    #= With Float32 floes, the device functions are compared with their Float64 results,
    which are the CPU functions' (see above). The CPU functions aren't a good reference for
    Float32 floes: they clip them partly in Float64 and return the overlap regions with
    Float64 points, so that, for example, `_many_intersect_normal_force!` finds none of a
    region's edges, while the device functions find most of them. In regions with two
    intersection points, the forces agree to about 1e-3; overlap areas of a few thousand m²
    lose precision in Float32 at these coordinates (up to about 65 km). In regions with more
    intersection points they don't agree, because `_many_intersect_normal_force` misses edges
    in Float32 (see its TODO). =#
    @testset "Pairs of floes with Float32" begin
        function _device_pairs(::Type{FT}) where FT
            floes = _collision_test_floes(FT)
            fwf = Subzero.FixedWidthFloes(floes)
            pairs = _test_pairs(length(floes))
            ws, out = _buffers(fwf, length(pairs))
            results = [_pair_run(out, ws, fwf, k, i, j, Constants{FT}(), _max_overlap(FT))
                for (k, (i, j)) in enumerate(pairs)]
            return floes, pairs, results, out
        end
        floes64, pairs, results64, out64 = _device_pairs(Float64)
        _, _, results32, out32 = _device_pairs(Float32)
        @test all(r -> r[1] >= 0, results32)
        # The same fuse decisions and numbers of forces
        @test results32 == results64
        # Relative differences of the forces, for regions with two and with more intersection
        # points
        Δ_two, Δ_many = Float64[], Float64[]
        for (k, (i, j)) in enumerate(pairs)
            results64[k][1] > 0 || continue
            regions = _regions_and_points(floes64, i, j)
            for r in 1:results64[k][1]
                f64, f32 = out64.force[k, r, :], Float64.(out32.force[k, r, :])
                Δ = iszero(f64) ? (iszero(f32) ? 0.0 : Inf) : hypot((f32 .- f64)...) / hypot(f64...)
                # the region of this force: the one with its area
                _, m = regions[argmin([abs(Subzero.area_poly(x, Float64) - out64.overlap[k, r])
                    for (x, _) in regions])]
                push!(m == 2 ? Δ_two : Δ_many, Δ)
            end
        end
        @test length(Δ_two) > 100
        @test count(<=(0.01), Δ_two) >= 0.95 * length(Δ_two)
        @test length(Δ_many) > 10
        @test_broken all(<=(0.01), Δ_many)
    end

    # The direction of the force from a region's edges, for regions with other than two
    # intersection points
    @testset "Force direction from many intersection points" begin
        FT = Float64
        floes = _collision_test_floes(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        # pairs with a region whose Δl or direction differs from the CPU function's
        n_tested, different = 0, Tuple{Int, Int}[]
        for (i, j) in _test_pairs(length(floes)), (region, m) in _regions_and_points(floes, i, j)
            _many_points(m) || continue
            cpu_dir = zeros(FT, 2)
            cpu_Δl = Subzero._many_intersect_normal_force!(FT, cpu_dir, region, floes.poly[i],
                FT(1e6))
            ring = [(GI.x(p), GI.y(p)) for p in GI.getpoint(GI.getexterior(region))][1:end-1]
            _, (dev_region,) = _array_polys([ring])
            Δl, dir = Subzero._many_intersect_normal_force(FT, dev_region,
                Subzero._floe_poly(fwf, i), FT(1e6))
            isequal((Δl, dir), (cpu_Δl, (cpu_dir[1], cpu_dir[2]))) || push!(different, (i, j))
            n_tested += 1
        end
        @test isempty(different)
        @test n_tested > 10
    end

    #= A narrow floe across a wider one: the edges of the overlap region on the first floe are
    on opposite sides, so their normals cancel out and the force has no direction (NaN). The
    CPU functions then throw an error in the clipping; the device functions return
    `GOG.CLIP_FAILED`. With `max_overlap = 1`, the floes don't fuse. =#
    @testset "No force direction with $FT" for FT in (Float64, Float32)
        square = [(0.0, 0.0), (1000.0, 0.0), (1000.0, 1000.0), (0.0, 1000.0)]
        bar = [(-200.0, 400.0), (1200.0, 400.0), (1200.0, 600.0), (-200.0, 600.0)]
        fs = FloeSettings(FT; min_floe_area = 0)
        floes = StructArray([Floe{FT}([[[FT(x), FT(y)] for (x, y) in [p; [p[1]]]]], FT(0.5);
            floe_settings = fs, rng = Xoshiro(1)) for p in (square, bar)])
        fwf = Subzero.FixedWidthFloes(floes)
        ws, out = _buffers(fwf, 1)
        consts = Constants{FT}()
        @test _pair_run(out, ws, fwf, 1, 1, 2, consts, FT(1)) == (GOG.CLIP_FAILED, false)
        @test_throws InexactError _cpu_pair_forces(floes, 1, 2, consts, FT(1))
    end

    # Buffers that are too small are reported, not overrun
    @testset "Too little room with $FT" for FT in (Float64, Float32)
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = _test_pairs(length(floes))
        ws, out = _buffers(fwf, length(pairs))
        counts = [first(_pair_run(out, ws, fwf, k, i, j, consts, max_overlap))
            for (k, (i, j)) in enumerate(pairs)]
        # the pair with the most forces
        i, j = pairs[argmax(counts)]
        n = maximum(counts)
        @test n > 1
        run_with(ws, out) = first(_pair_run(out, ws, fwf, 1, i, j, consts, max_overlap))
        big_ws, big_out = _buffers(fwf, 1)
        @test run_with(big_ws, big_out) == n
        @test run_with(big_ws, Subzero.PairForces(FT, 1, n - 1)) == GOG.CLIP_TOO_SMALL
        @test run_with(first(_buffers(fwf, 1; max_ipoints = 1)), big_out) == GOG.CLIP_TOO_SMALL
        @test run_with(first(_buffers(fwf, 1; max_regions = n - 1)), big_out) ==
            GOG.CLIP_TOO_SMALL
    end

    # Kernels: one thread per pair, against the same function called on the CPU without a kernel
    function _pair_forces_kernel!(k, floes, pairs, ws, out, results, consts, Δt, max_overlap)
        n, fuse = Subzero.floe_floe_interaction!(Subzero.PairForces(out, k),
            Subzero.CollisionWorkspace(ws, k), floes, pairs[k, 1], pairs[k, 2], consts, Δt,
            max_overlap)
        results[k, 1] = n
        results[k, 2] = fuse
        return
    end
    #= With Float32 floes, the kernel doesn't compute in Float64. Checked without compiling
    it with POCL, which takes minutes (see `_kernel_computes_in_float64`). =#
    @testset "No Float64 maths in the kernel with Float32" begin
        FT = Float32
        floes = _collision_test_floes(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = Int32.(reduce(vcat, [[i j] for (i, j) in _test_pairs(length(floes))]))
        npairs = size(pairs, 1)
        args = (fwf, pairs, _buffers(fwf, npairs)..., zeros(Int32, npairs, 2), Constants{FT}(),
            Δt, _max_overlap(FT))
        @test !_kernel_computes_in_float64(() ->
            _gog_launch!(_pair_forces_kernel!, Subzero.KernelAbstractions.CPU(), npairs, args...))
    end

    #= Float64 on the CPU backend and the GPUs, Float32 on the GPUs only: POCL takes minutes
    to compile this kernel. JLBackend compiles with POCL too, without reusing CPU()'s compile,
    so it doesn't run it. With Float32, the kernel's maths is checked above, without
    compiling it. =#
    _kernel_float_types(backend) =
        backend isa JLArrays.JLBackend ? () :
        backend isa typeof(Subzero.KernelAbstractions.CPU()) ? (Float64,) : (Float64, Float32)
    @testset "Kernels on $backend with $FT" for backend in test_backends(),
            FT in _kernel_float_types(backend)
        floes = _collision_test_floes(FT)
        consts = Constants{FT}()
        max_overlap = _max_overlap(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        pairs = Int32.(reduce(vcat, [[i j] for (i, j) in _test_pairs(length(floes))]))
        npairs = size(pairs, 1)
        outputs() = (_buffers(fwf, npairs)..., zeros(Int32, npairs, 2))
        ref = outputs()
        for k in 1:npairs
            _pair_forces_kernel!(k, fwf, pairs, ref..., consts, Δt, max_overlap)
        end
        @test all(>=(0), ref[3][:, 1])
        dev = Subzero.adapt(backend, outputs())
        _gog_launch!(_pair_forces_kernel!, backend, npairs, Subzero.adapt(backend, fwf),
            Subzero.adapt(backend, pairs), dev..., consts, Δt, max_overlap)
        got = Subzero.adapt(Array, dev)
        # pairs whose results differ from those on the CPU without a kernel
        different = [Tuple(pairs[k, :]) for k in 1:npairs if got[3][k, :] != ref[3][k, :] ||
            !isequal(_device_rows(got[2], k, ref[3][k, 1]), _device_rows(ref[2], k, ref[3][k, 1]))]
        @test isempty(different)
        @test sum(ref[3][:, 1]) > 50
    end
end
