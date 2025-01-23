using JLD2, Logging, NCDatasets, Random, SplitApplyCombine,
    Statistics, StructArrays, Subzero, VoronoiCells
import GeometryOps as GO
import GeometryOps.GeoInterface as GI
using Test

@testset "Subzero.jl" begin
    include("test_simulation_components/grids.jl")
    include("test_simulation_components/domain_components/boundaries.jl")
    include("test_simulation_components/domain_components/domains.jl")
    include("test_simulation_components/domain_components/topography.jl")
    include("test_simulation_components/test_stress_calculators.jl")
    include("test_simulation_components/oceans.jl")
    include("test_simulation_components/atmos.jl")
    include("test_physical_processes/test_update_floe.jl")
    include("test_physical_processes/test_collisions.jl")
    include("test_physical_processes/test_coupling.jl")
    include("test_physical_processes/test_process_settings.jl")
    include("test_physical_processes/test_fractures.jl")
    include("test_physical_processes/test_simplification.jl")
    include("test_physical_processes/test_ridge_raft.jl")
    include("test_physical_processes/test_welding.jl")
    include("test_floe.jl")
    include("test_floe_utils.jl")
    include("test_model.jl")
    include("test_output.jl")
    include("test_simulation.jl")
    include("test_conservation.jl")
end
