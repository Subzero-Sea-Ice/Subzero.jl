#=
Floes for testing `add_ghosts!` on the device against the CPU version, in a 10 km × 10 km
domain: random star-shaped floes in and around the domain, so that many cross a boundary or
a corner and some have their centroid outside of the domain; rectangles that only touch a
boundary, cross a corner, or have their centroid outside of the domain; and inactive floes
that cross a boundary.
=#
function _ghost_test_floes(::Type{FT}; n = 150) where FT
    rng = Xoshiro(5)
    floe_settings = FloeSettings(FT; min_floe_area = 0)
    floes = StructArray{Floe{FT}}(undef, 0)
    function add!(points)
        ring = [[Float64(x), Float64(y)] for (x, y) in points]
        push!(ring, ring[1])
        push!(floes, Floe{FT}([ring], 0.5 + rand(rng); floe_settings, rng))
    end
    rectangle(x1, x2, y1, y2) = [(x1, y1), (x2, y1), (x2, y2), (x1, y2)]
    add!(rectangle(9e3, 1e4, 2e3, 3e3))  # touches the east boundary
    add!(rectangle(0, 1e3, 5e3, 6e3))  # touches the west boundary
    add!(rectangle(4e3, 5e3, 9e3, 1e4))  # touches the north boundary
    add!(rectangle(4e3, 5e3, 0, 1e3))  # touches the south boundary
    add!(rectangle(9e3, 1e4, 9e3, 1e4))  # touches the north-east corner
    add!(rectangle(9.5e3, 1.05e4, 9.5e3, 1.05e4))  # crosses the north-east corner
    add!(rectangle(-5e2, 5e2, -5e2, 5e2))  # crosses the south-west corner
    add!(rectangle(9.8e3, 1.1e4, 2e3, 3e3))  # centroid east of the domain
    add!(rectangle(4e3, 5e3, -1.2e3, 2e2))  # centroid south of the domain
    add!(rectangle(-1.2e3, 2e2, 9.6e3, 1.1e4))  # centroid north-west of the domain
    for _ in 1:n
        cx, cy = -1e3 .+ 1.2e4 .* rand(rng, 2)
        r = 300 + 1200rand(rng)
        θs = sort(2π .* rand(rng, rand(rng, 5:10)))
        add!([(cx + r * (0.5 + rand(rng) / 2) * cos(θ), cy + r * (0.5 + rand(rng) / 2) * sin(θ)) for θ in θs])
    end
    # inactive floes that cross a boundary
    floes.status[6].tag = Subzero.remove
    floes.status[8].tag = Subzero.fuse
    return floes
end

#=
Add ghosts to floes on the CPU and on `backend`, and return the CPU result and the
device result (as a FixedWidthFloes on the CPU).
=#
function _ghosts_on(backend, floes, domain; capacity = length(floes) + cld(length(floes), 4), workspace = nothing)
    cpu = deepcopy(floes)
    Subzero.add_ghosts!(cpu, domain)
    dev = Subzero.adapt(backend, Subzero.FixedWidthFloes(floes; capacity))
    dev = Subzero.add_ghosts!(dev, domain, backend; workspace)
    return cpu, Subzero.adapt(Array, dev)
end

# Remove the ghosts from floes on the CPU, as timestep_sim! does
function _remove_ghosts_cpu!(floes, n)
    for i in reverse((n + 1):length(floes))
        StructArrays.foreachfield(field -> deleteat!(field, i), floes)
    end
    empty!.(floes.ghosts)
    return floes
end

# The fields of the FixedWidthFloes `got` that differ from those of the floes `expected`
function _differing_fields(got, expected::StructArray)
    want = Subzero.FixedWidthFloes(expected)
    n = length(expected)
    got.n_floes == n || return [:n_floes]
    return filter(fieldnames(Subzero.FixedWidthFloes)) do field
        field == :n_floes && return false
        a, b = getfield(got, field), getfield(want, field)
        size(a)[2:end] != size(b)[2:end] || selectdim(a, 1, 1:n) != selectdim(b, 1, 1:n)
    end
end

@testset "Ghosts on the device" begin
    @testset "add_ghosts! on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        grid = RegRectilinearGrid(FT; x0 = 0.0, xf = 1e4, y0 = 0.0, yf = 1e4, Δx = 1e3, Δy = 1e3)
        periodic(D) = PeriodicBoundary(D, FT; grid)
        collision(D) = CollisionBoundary(D, FT; grid)
        domains = (
            both = Domain(; north = periodic(North), south = periodic(South),
                east = periodic(East), west = periodic(West)),
            east_west = Domain(; north = collision(North), south = collision(South),
                east = periodic(East), west = periodic(West)),
            north_south = Domain(; north = periodic(North), south = periodic(South),
                east = collision(East), west = collision(West)),
            none = Domain(; north = collision(North), south = collision(South),
                east = collision(East), west = collision(West)),
        )
        floes = _ghost_test_floes(FT)
        n = length(floes)
        for (name, domain) in pairs(domains)
            # with the default capacity, and without spare capacity, so that the arrays have to grow
            for capacity in (n + cld(n, 4), n)
                cpu, got = _ghosts_on(backend, floes, domain; capacity)
                @test isempty(_differing_fields(got, cpu))
                if name == :both
                    # the test floes have floes with 0, 1, and 3 ghosts, and inactive ones
                    @test Set(length.(cpu.ghosts[1:n])) == Set([0, 1, 3])
                    @test isempty(cpu.ghosts[6]) && isempty(cpu.ghosts[8])
                    @test length(cpu) > 1.25n  # more ghosts than spare rows
                elseif name == :none
                    @test got.n_floes == n
                end
            end
        end
        # A workspace can be reused, and the arrays are shared when the ghosts fit
        fwf = Subzero.adapt(backend, Subzero.FixedWidthFloes(floes; capacity = 2n))
        workspace = Subzero.GhostWorkspace(fwf)
        for _ in 1:2
            dev = Subzero.add_ghosts!(Subzero.set_n_floes(fwf, n), domains.east_west, backend; workspace)
            @test dev.poly === fwf.poly
            cpu = deepcopy(floes)
            Subzero.add_ghosts!(cpu, domains.east_west)
            @test isempty(_differing_fields(Subzero.adapt(Array, dev), cpu))
            # remove the ghosts again, as timestep_sim! does
            fwf = Subzero.adapt(backend, Subzero.FixedWidthFloes(floes; capacity = 2n))
        end
        #= Removing the ghosts, as timestep_sim! does, and adding them again to the same
        arrays, as the next timestep would, when the rows of the old ghosts still hold their
        values. =#
        for (name, domain) in pairs(domains)
            cpu = deepcopy(floes)
            Subzero.add_ghosts!(cpu, domain)
            _remove_ghosts_cpu!(cpu, n)
            dev = Subzero.adapt(backend, Subzero.FixedWidthFloes(floes))
            dev = Subzero.add_ghosts!(dev, domain, backend)
            dev = Subzero.remove_ghosts!(dev, n, backend)
            Subzero.KernelAbstractions.synchronize(backend)
            @test isempty(_differing_fields(Subzero.adapt(Array, dev), cpu))
            # some floes swapped places with a ghost and keep their new place
            name == :both && @test cpu.centroid != floes.centroid
            Subzero.add_ghosts!(cpu, domain)
            dev = Subzero.add_ghosts!(dev, domain, backend)
            @test isempty(_differing_fields(Subzero.adapt(Array, dev), cpu))
        end
        # No floes
        empty_floes = Subzero.adapt(backend, Subzero.FixedWidthFloes(floes[1:0]))
        @test Subzero.add_ghosts!(empty_floes, domains.both, backend).n_floes == 0
    end
end
