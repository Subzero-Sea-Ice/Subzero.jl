@testset "Collisions" begin
    FT = Float64
    Δt = 10
    Δh = 0.0
    consts = Constants()
    collision_settings = CollisionSettings()
    spinlock = Threads.SpinLock()
    supress_warnings = true
    # Grid setup
    Lx = 1e5
    Ly = Lx
    grid = RegRectilinearGrid(; x0 = -Lx, xf = Lx, y0 = -Ly, yf = Ly, Δx = 1e4, Δy = 1e4)

    # Boundary setup
    p_n_bound = PeriodicBoundary(North; grid)
    p_s_bound = PeriodicBoundary(South; grid)
    p_e_bound = PeriodicBoundary(East; grid)
    p_w_bound = PeriodicBoundary(West; grid)

    c_n_bound =  CollisionBoundary(North; grid)
    c_s_bound = CollisionBoundary(South; grid)
    c_e_bound = CollisionBoundary(East; grid)
    c_w_bound = CollisionBoundary(West; grid)

    o_n_bound = OpenBoundary(North; grid)
    o_s_bound = OpenBoundary(South; grid)
    o_e_bound = OpenBoundary(East; grid)
    o_w_bound = OpenBoundary(West; grid)

    topos = initialize_topography_field(FT; coords = [[[[1e4, 0.0], [0.0, 1e4], [1e4, 2e4], [2e4, 1e4], [1e4, 0.0]]]])
        
    topo_domain = Domain(; north = p_n_bound, south = p_s_bound, east = c_e_bound, west = o_w_bound, topography = topos)
    collision_domain = Domain(; north = c_n_bound, south = c_s_bound, east = c_e_bound, west = c_w_bound)
    open_domain = Domain(; north = o_n_bound, south = o_s_bound, east = o_e_bound, west = o_w_bound)
    ew_periodic_domain = Domain(; north = o_n_bound, south = o_s_bound, east = p_e_bound, west = p_w_bound)
    ns_periodic_domain = Domain(; north = p_n_bound, south = p_s_bound, east = o_e_bound, west = o_w_bound)
    double_periodic_domain = Domain(; north = p_n_bound, south = p_s_bound, east = p_e_bound, west = p_w_bound)

    @testset "Floe-Floe Interactions" begin
        hmean = 0.25
        max_overlap = 0.55
        # Floe coordinates
        tri_coord = [[[0.0, 0.0], [1e4, 3e4], [2e4, 0], [0.0, 0.0]]]
        corner_rect_coord = [[[0.0, 2.5e4], [0.0, 2.9e4], [2e4, 2.9e4], [2e4, 2.5e4], [0.0, 2.5e4]]]
        small_shift_corner_rect_coord = Subzero.translate(corner_rect_coord, 0.5e4, 0.0)
        big_shift_corner_rect_coord = Subzero.translate(corner_rect_coord, 1.9999999e4, 0.0)
        middle_rect_coord = [[ [1.8e4, 2.7e4], [1.8e4, 2.8e4], [2.1e4, 2.8e4], [2.1e4, 2.7e4], [1.8e4, 2.7e4]]]
        cshape_coord = [[[0.5e4, 2.7e4], [0.5e4, 3.5e4], [1.5e4, 3.5e4], [1.5e4, 2.7e4], [1.25e4, 2.7e4], [1.25e4, 3e4], [1e4, 3e4], [1e4, 2.7e4], [0.5e4, 2.7e4]]]

        # Test interactions between triange tip intersecting with a rectangle
        tri, corner_rect = Floe(tri_coord, hmean, Δh), Floe(corner_rect_coord, hmean, Δh)
        tri.u, corner_rect.v = 0.1, -0.1
        Subzero.floe_floe_interaction!(tri, 1, corner_rect, 2, consts, Δt, max_overlap)
        @test isapprox(tri.interactions[1, xforce], -64613382.47, atol = 1e-2)
        @test isapprox(tri.interactions[1, yforce], -521498991.51, atol = 1e-2)
        @test isapprox(tri.interactions[1, xpoint], 10000.00, atol = 1e-2)
        @test isapprox(tri.interactions[1, ypoint], 26555.55, atol = 1e-2)
        @test isapprox(tri.interactions[1, overlap], 8000000, atol = 1e-2)
        @test tri.status.tag != Subzero.fuse && corner_rect.status.tag != Subzero.fuse
        @test isempty(tri.status.fuse_idx)
        Subzero.calc_torque!(tri)
        @test isapprox(tri.interactions[1, torque], 1069710443203.99, atol = 1e-2)

        # Sideways C intersected with rectangle so there are two areas of overlap
        cshape_floe, corner_rect = Floe(cshape_coord, hmean, Δh), Floe(corner_rect_coord, hmean, Δh)
        cshape_floe.u, corner_rect.v = 0.3, -0.1
        Subzero.floe_floe_interaction!(cshape_floe, 1, corner_rect, 2, consts, Δt, max_overlap)
        @test isapprox(cshape_floe.interactions[1, xforce], -163013665.41, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, xforce], -81506832.70, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, yforce], 804819565.60, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, yforce], 402409782.80, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, xpoint], 7500.00, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, xpoint], 13750.00, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, ypoint], 28000.00, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, ypoint], 28000.00, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, overlap], 10000000, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, overlap], 5000000, atol = 1e-2)
        @test tri.status.tag != Subzero.fuse && corner_rect.status.tag != Subzero.fuse
        Subzero.calc_torque!(cshape_floe)
        @test isapprox(cshape_floe.interactions[1, torque], -2439177121266.03, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, torque], 1295472581868.05, atol = 1e-2)

        # Floes overlapping more than 55%  - rectangle and shifted rectangle
        corner_rect, small_shift_corner_rect = Floe(corner_rect_coord, hmean, Δh), Floe(small_shift_corner_rect_coord, hmean, Δh)
        corner_rect.v, small_shift_corner_rect.v = -0.1, -0.1
        Subzero.floe_floe_interaction!(corner_rect, 1, small_shift_corner_rect, 2, consts, Δt, max_overlap)
        @test corner_rect.status.tag == Subzero.fuse 
        @test corner_rect.status.fuse_idx == [2]
        @test isempty(corner_rect.interactions)

        # Two floes (original and original with a small shift) overlapping by more than 55%
        corner_rect, middle_rect = Floe(corner_rect_coord, hmean, Δh), Floe(middle_rect_coord, hmean, Δh)
        corner_rect.v = -0.1
        Subzero.floe_floe_interaction!(corner_rect, 1, middle_rect, 2, consts, Δt, max_overlap)
        @test corner_rect.status.tag == Subzero.fuse 
        @test corner_rect.status.fuse_idx == [2]
        
        # Two floes (original and original with a big shift) have such a small overlap that forces are not calculated
        corner_rect, big_shift_corner_rect = Floe(corner_rect_coord, hmean, Δh), Floe(big_shift_corner_rect_coord, hmean, Δh)
        corner_rect.v, big_shift_corner_rect.v = -0.1, -0.1
        Subzero.floe_floe_interaction!(big_shift_corner_rect, 1, corner_rect, 2, consts, Δt, max_overlap)
        @test isempty(big_shift_corner_rect.interactions)
    end

    @testset "Floe Boundary Interactions" begin
        hmean = 0.25
        max_overlap = 0.75
        ## Floe coordinates
        # Floe overlapping with northern periodic boundary
        north_coords = [[[5e4, 9.75e4], [5e4, 10.05e4], [7e4, 10.05e4], [7e4, 9.75e4], [5e4, 9.75e4]]]
        # Diagonal floe barely overlaping with eastern collision boundary
        east_small_coords = [[[9.5e4, 0.0], [9e4, 0.5e4], [10e4, 2.5e4], [10.05e4, 2e4], [9.5e4, 0.0]]]
        # Floe overlapping with eastern collision boundary by more than 75% to trigger overlap condition
        east_large_coords = [[[9e4, -7e4], [9e4, -5e4], [1.4e5, -5e4], [1.4e5, -7e4], [9e4, -7e4]]]
        # Floe overlapping with western open boundary
        west_coords = [[[-9.75e4, 7e4], [-9.75e4, 5e4], [-10.05e4, 5e4], [-10.05e4, 7e4], [-9.75e4, 7e4]]]
        # Floe overlapping with boundary with more than one region
        cshape_coord = [[[9.5e4, 7e4], [9.5e4, 9e4], [1.05e5, 9e4], [1.05e5, 8.5e4], [9.9e4, 8.5e4], [9.9e4, 8e4], [1.05e5, 8e4], [1.05e5, 7e4], [9.5e4, 7e4]]]
        # Floe overlapping with topography element
        topo_overlap_coords = [[[-0.5e4, 0.0], [-0.5e4, 0.75e4], [0.5e4, 0.75e4], [0.5e4, 0.0], [-0.5e4, 0.0]]]
        # Floe in a corner hitting more than one wall at once
        corner_coords = [[[9.5e4, 7e4], [9e4, 7.5e4], [10e4, 1.05e5], [10.05e4, 9.5e4], [9.5e4, 7e4]]]
    
        # Test floe overlapping slightly with collision boundary (one region)
        east_small_floe = Floe(east_small_coords, hmean, Δh)
        east_small_floe.u, east_small_floe.v = 0.5, 0.25
        Subzero.floe_domain_interaction!(east_small_floe, topo_domain, consts, Δt, max_overlap)
        @test east_small_floe.interactions[1, floeidx] == -3
        @test isapprox(east_small_floe.interactions[1, xforce], -311304795.629, atol = 1e-3)
        @test isapprox(east_small_floe.interactions[1, yforce], -23618874.648, atol = 1e-3)
        @test isapprox(east_small_floe.interactions[1, overlap], 1704545.454, atol = 1e-3)
        @test isapprox(east_small_floe.interactions[1, xpoint], 100166.666, atol = 1e-3)
        @test isapprox(east_small_floe.interactions[1, ypoint], 21060.606, atol = 1e-3)
        
        # Test floe overlapping slightly with collision boundary (two regions)
        cshape_floe = Floe(cshape_coord, hmean, Δh)
        cshape_floe.v = -0.1
        Subzero.floe_domain_interaction!(cshape_floe, topo_domain, consts, Δt, max_overlap)
        @test cshape_floe.interactions[1, floeidx] == -3
        @test cshape_floe.interactions[2, floeidx] == -3
        @test isapprox(cshape_floe.interactions[1, xforce], -2876118708.17, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, xforce], -5752237416.35, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, yforce], 575223741.63, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, yforce], 1150447483.27, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, xpoint], 102500, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, xpoint], 102500, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, ypoint], 87500, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, ypoint], 75000, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[1, overlap], 25000000, atol = 1e-2)
        @test isapprox(cshape_floe.interactions[2, overlap], 50000000, atol = 1e-2)

        # Test floe overlapping >75% with collision boundary
        east_large_floe = Floe(east_large_coords, hmean, Δh)
        east_large_floe.u, east_large_floe.v = -0.4, 0.2
        Subzero.floe_domain_interaction!(east_large_floe, topo_domain, consts, Δt, max_overlap)
        @test isempty(east_large_floe.interactions)
        @test east_large_floe.status.tag == Subzero.remove

        # Test floe overlapping >75% with collision boundary -> but max overlap is now 1
        east_large_floe = Floe(east_large_coords, hmean, Δh)
        east_large_floe.u, east_large_floe.v = -0.4, 0.2
        Subzero.floe_domain_interaction!(east_large_floe, topo_domain, consts, Δt, 1.0)
        @test !isempty(east_large_floe.interactions)
        @test east_large_floe.num_inters > 0

        # Test floe passing through open boundary is killed
        west_floe = Floe(west_coords, hmean, Δh)
        Subzero.floe_domain_interaction!(west_floe, topo_domain, consts, Δt, max_overlap)
        @test west_floe.status.tag == Subzero.remove

        # Test floes not not interact with periodic boundary
        north_floe = Floe(north_coords, hmean, Δh)
        Subzero.floe_domain_interaction!(north_floe, topo_domain, consts, Δt, max_overlap)
        @test north_floe.status.tag == north_floe.status.tag && north_floe.interactions == north_floe.interactions

        # Test floe overlapping with topography -> different from Subzero since topography are now boundaries
        topo_overlap_floe = Floe(topo_overlap_coords, hmean, Δh)
        Subzero.floe_domain_interaction!(topo_overlap_floe, topo_domain, consts, Δt, max_overlap)
        @test topo_overlap_floe.interactions[1, xforce] < 0
        @test topo_overlap_floe.interactions[1, yforce] < 0
        @test topo_overlap_floe.interactions[1, floeidx] == -5

        # Test floe in corner hitting more than one wall at once
        corner_floe = Floe(corner_coords, hmean, Δh)
        Subzero.floe_domain_interaction!(corner_floe, collision_domain, consts, Δt, max_overlap)
        @test all(corner_floe.interactions[:, xforce] .<= 0)
        @test all(corner_floe.interactions[:, yforce] .<= 0)
    end
    
    @testset "Add Ghosts" begin
        hmean = 0.5
        ## Floe coordinates
        # corner floe - overlaps with north and east boundary
        coords1 = [[[9.9e4, 9.9e4], [9.9e4, 1.02e5], [1.02e5, 1.02e5], [1.02e5, 9.9e4], [9.9e4, 9.9e4]]]
        # overlaps with western boundary
        coords2 = [[[-1.01e5, 7e4], [-1.01e5, 8e4], [-8e4, 8e4], [-8e4, 7e4], [-1.01e5, 7e4]]]
        # overlaps with northern boundary
        coords3 = [[[-2e4, 9.5e4], [-2e4, 1.1e5], [-1e4, 1.1e5], [-1e4, 9.5e4], [-2e4, 9.5e4]]]
        # doesn't overlap with any boundary
        coords4 = [[[0.0, 0.0], [0.0, 2e4], [2e4, 2e4], [2e4, 0.0], [0.0, 0.0]]]
        # List of above coordinates
        coord_list = [coords1, coords2, coords3, coords4]

        # Make sure nothing is added with non-periodic domain
        floe_arr = initialize_floe_field(FT, coord_list, open_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, open_domain)
        @test floe_arr.coords == coord_list

        # Add ghost floes in east-west direction
        floe_arr = initialize_floe_field(FT, coord_list, ew_periodic_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, ew_periodic_domain)
        @test -1e5 < floe_arr[1].centroid[1] < 1e5
        @test -1e5 < floe_arr[2].centroid[2] < 1e5
        @test floe_arr.coords[1] == Subzero.translate(coords1, -2e5, 0.0)
        @test floe_arr.coords[2:4] == coord_list[2:4]
        @test floe_arr.coords[5] == coords1
        @test floe_arr.coords[6] == Subzero.translate(coords2, 2e5, 0.0)
        @test floe_arr.id == [1, 2, 3, 4, 1, 2]
        @test floe_arr.ghost_id == [0, 0, 0, 0, 1, 1]
        @test floe_arr.ghosts[1] == [5]
        @test floe_arr.ghosts[2] == [6]
        @test floe_arr.ghosts[3:6] == [[], [], [], []]

        # Add ghost floes in the north-south direction
        floe_arr = initialize_floe_field(FT, coord_list, ns_periodic_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, ns_periodic_domain)
        @test -1e5 < floe_arr[1].centroid[2] < 1e5
        @test -1e5 < floe_arr[3].centroid[2] < 1e5
        @test floe_arr.coords[1] == Subzero.translate(coords1, 0.0, -2e5)
        @test floe_arr.coords[3] == Subzero.translate(coords3, 0.0, -2e5)
        @test floe_arr.coords[[2, 4]] == coord_list[[2, 4]]
        @test floe_arr.coords[5] == coords1
        @test floe_arr.coords[6] == coords3
        @test floe_arr.id == [1, 2, 3, 4, 1, 3]
        @test floe_arr.ghost_id == [0, 0, 0, 0, 1, 1]
        @test floe_arr.ghosts[1] == [5]
        @test floe_arr.ghosts[3] == [6]
        @test floe_arr.ghosts[[2; 4:6]] == [[], [], [], []]

        # Add ghosts in both east-west and north-south directions
        floe_arr = initialize_floe_field(FT, coord_list, double_periodic_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, double_periodic_domain)
        @test -1e5 < floe_arr.centroid[1][1] < 1e5
        @test -1e5 < floe_arr.centroid[1][2] < 1e5
        @test floe_arr.coords[1] == Subzero.translate(coords1, -2e5, -2e5)
        @test floe_arr.coords[3] == Subzero.translate(coords3, 0.0, -2e5)
        @test floe_arr.coords[[2, 4]] == coord_list[[2, 4]]
        @test floe_arr.coords[5] == coords1
        @test floe_arr.coords[6] == Subzero.translate(coords2, 2e5, 0.0)
        @test floe_arr.coords[7] == Subzero.translate(coords1, 0.0, -2e5)
        @test floe_arr.coords[8] == Subzero.translate(coords1, -2e5, 0.0)
        @test floe_arr.coords[9] == coords3
        @test floe_arr.id == [1, 2, 3, 4, 1, 2, 1, 1, 3]
        @test floe_arr.ghost_id == [0, 0, 0, 0, 1, 1, 2, 3, 1]
        @test floe_arr.ghosts[1] == [5, 7, 8]
        @test floe_arr.ghosts[2] == [6]
        @test floe_arr.ghosts[3] == [9]
        @test floe_arr.ghosts[4:9] == [[], [], [], [], [], []]
    end
    @testset "Ghost Collisions" begin
        hmean = 0.5
        ## Floe coordinates
        # L-shaped polygon extending beyond top and right domain edges
        lshape_coords = [splitdims([Lx/2 Lx/2 3*Lx/4 3*Lx/4 Lx+10000 Lx+10000; Ly/2 Ly+10000 Ly+10000 3*Ly/4 3*Ly/4 Ly/2])]
        # Oval floe touching two edges of lshape_coords
        th, r = 0:pi/50:2*pi, Ly/4+1000
        oval_coords = [invert([r * cos.(th) .+ (Lx-1), r * sin.(th) .+ (Ly-1)])]
        # Tall rectangle in top right corner of domain
        tall_rect_coords = [splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))]
        # Long rectangle in bottom left corner of domain (ghost hits tall rectangle)
        long_rect_coords = [splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))]
        # Shifted tall rectangle to ghost position in bottom right corner
        shifted_down_tall_rect_coords = Subzero.translate(tall_rect_coords, 0.0, -2Ly)
        # Shifted tall rectangle to ghost postion in top left corner
        shifted_left_tall_rect_coords = Subzero.translate(tall_rect_coords, -2Lx, 0.0)
        # Shifted long rectangle to ghost postion in bottom right corner
        shifted_right_long_rect_coords = Subzero.translate(long_rect_coords, 2Lx, 0.0)
        # Shifted long rectangle to ghost postion in top left corner
        shifted_up_long_rect_coords = Subzero.translate(long_rect_coords, 0.0, 1.615Ly)
        # Small rectangle in the corner that has 3 ghosts in all other corners
        small_corner_rect_coords = [[[-1.1e5, -1.1e5], [-1.1e5, -9.5e4], [-9.5e4, -9.5e4], [-9.5e4, -1.1e5], [-1.1e5, -1.1e5]]]
        # triangle in the middle of the domain with no ghosts - touches 3/4 corners
        large_tri_coords = [[[-1e5, -1e5], [-1e5, 1e5], [1e5, -1e5], [-1e5, -1e5]]]
        # rectangle along south boundary, ghost along north boundary
        south_bound_rect_coords = [[[-9.8e4, -1.1e5], [-9.8e4, -9.5e4], [9.8e4, -9.5e4], [9.8e4, -1.1e5], [-9.8e4, -1.1e5]]]

        # Parent-parent collison (parents are touching) 
        floe_arr = initialize_floe_field(FT, [lshape_coords, oval_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        xforce_vals, yforce_vals = abs(floe_arr[1].collision_force[1]), abs(floe_arr[1].collision_force[2])
        f1_torque, f2_torque = floe_arr[1].collision_trq, floe_arr[2].collision_trq

        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        # 1 and 2 are the "parent" floes - floe 1 and floe 2 interact
        @test xforce_vals == abs(floe_arr[1].collision_force[1]) == abs(floe_arr[2].collision_force[1])
        @test yforce_vals == abs(floe_arr[2].collision_force[2]) == abs(floe_arr[2].collision_force[2])
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        # All other ghost floes aren't calculated
        @test first.(floe_arr[3:end].collision_force) == zeros(6)
        @test last.(floe_arr[3:end].collision_force) == zeros(6)
        @test floe_arr[3:end].collision_trq == zeros(6)

        # Ghost-Ghost collision (parents aren't touching, only ghosts touch)
        tall_rect_coords = [splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))]
        long_rect_coords = [splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))]
        floe_arr = initialize_floe_field(FT, [tall_rect_coords, long_rect_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        trans_arr = initialize_floe_field(FT, [shifted_down_tall_rect_coords, shifted_right_long_rect_coords], double_periodic_domain, hmean, Δh; supress_warnings)

        Subzero.timestep_collisions!(trans_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        xforce_vals = abs(trans_arr[1].collision_force[1])
        yforce_vals = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!( floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        # floes 1 and 2 are the parents - floe 4 is floe 1's ghost and floe 3 is floe 2's ghost - floe 3 and 4 collide 
        @test repeat([xforce_vals], 2) == abs.(first.(floe_arr.collision_force[1:2]))
        @test repeat([yforce_vals], 2) == abs.(last.(floe_arr.collision_force[1:2]))
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        # interactions copied from ghosts
        @test floe_arr[1].interactions[:, [1:5; 7]] == floe_arr[4].interactions[:, [1:5; 7]]
        @test floe_arr[2].interactions[:, [1:5; 7]] == floe_arr[3].interactions[:, [1:5; 7]]

        # Parent-Ghost Collision
        floe_arr = initialize_floe_field(FT, [tall_rect_coords, shifted_up_long_rect_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        trans_arr = initialize_floe_field(FT, [shifted_left_tall_rect_coords, shifted_up_long_rect_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        Subzero.timestep_collisions!(trans_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        xforce_vals = abs(trans_arr[1].collision_force[1])
        yforce_vals = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        # Floe 1's ghost if floe 4 and floe 2's ghost is floe 3 and floe 3 and floe 1 interact
        @test repeat([xforce_vals], 2) == abs.(first.(floe_arr.collision_force[1:2]))
        @test repeat([yforce_vals], 2) == abs.(last.(floe_arr.collision_force[1:2]))
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        @test floe_arr[2].interactions[:, [1:5; 7]] == floe_arr[3].interactions[:, [1:5; 7]]
        @test isempty(floe_arr[4].interactions)

        # Parent and ghosts hitting the same floe
        floe_arr =  initialize_floe_field(FT, [small_corner_rect_coords, large_tri_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, double_periodic_domain)
        @test length(floe_arr) == 5
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        @test size(floe_arr[1].interactions)[1] == 3
        @test size(floe_arr[2].interactions)[1] == 3
        @test floe_arr[1].interactions[1, Subzero.xforce] != floe_arr[1].interactions[2, Subzero.xforce] && floe_arr[1].interactions[1, Subzero.xforce] != floe_arr[1].interactions[3,Subzero.xforce]
        @test floe_arr[1].interactions[1, Subzero.yforce] != floe_arr[1].interactions[2, Subzero.yforce] && floe_arr[1].interactions[1, Subzero.yforce] != floe_arr[1].interactions[3, Subzero.yforce]

        floe_arr =  initialize_floe_field(FT, [small_corner_rect_coords, south_bound_rect_coords], double_periodic_domain, hmean, Δh; supress_warnings)
        add_ghosts!(floe_arr, double_periodic_domain)
        @test length(floe_arr) == 6
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, consts, Δt, collision_settings, spinlock)
        @test size(floe_arr[1].interactions)[1] == 2
        @test size(floe_arr[2].interactions)[1] == 2
        @test floe_arr[1].interactions[1, Subzero.xpoint] != floe_arr[1].interactions[2, Subzero.xpoint]
        @test floe_arr[1].interactions[1, Subzero.ypoint] == floe_arr[1].interactions[2, Subzero.ypoint]
    end
end