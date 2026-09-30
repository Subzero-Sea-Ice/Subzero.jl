# # Ice-ocean Coupled Simulation with Oceananigans

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../coupled_ice_ocean/oceananigans_coupled_example.mp4" type="video/mp4">
# </video>
# ```

# This simulation provides an example of mechanically coupling Subzero.jl to the ocean large-eddy  
# simulation model, [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/). 
# The ocean model is set up to simulate the Eady instability problem in a doubly-periodic domain. 
# As the ocean current fields evolve, they provide forcing to Subzero.jl. In turn, the 
# ice-ocean stresses computed by Subzero.jl are passed back to Oceananigans as surface stress 
# boundary conditions.

# !!! note
#       Note that this example uses CPU execution for simplicity, but Oceananigans.jl is capable of running on GPUs.
#       This is still possible in the coupled framework, but requires a few additional setup steps.

using Subzero
using Oceananigans, Oceananigans.Units
using CairoMakie, GeoInterfaceMakie
using JLD2, Random, Printf

# ## User Inputs

# Grid constants:
const FT = Float64          # Float type used to run simulation
const Lx = 51200            # grid x-length
const Ly = 51200            # grid y-length
const Lz = 16               # grid z-length
const Δgrid = 400           # grid cell edge-size
const Δz = 2meter           # grid cell vertical edge-size
const Nx = Int(Lx/Δgrid)    # number of grid cells in x-direction
const Ny = Int(Ly/Δgrid)    # number of grid cells in y-direction
const Nz = Int(Lz/Δz)       # number of grid cells in z-direction

# Coupling between the two models can be done at the ocean timestep (as done here), or at some 
# multiple (an approximation that reduces computational cost at the expense of accuracy and stability).
# The sea ice model often needs a much faster timestep than the ocean model in order to resolve collisions.
const simTime = 10days              # total simulation time
const Δtᵢ = 5seconds                # timestep for sea ice
const Δtₒ = 150seconds              # timestep for ocean
const ΔtCpl = Δtₒ                   # timestep for ice-ocean coupling
const nΔtᵢ = Int(simTime/Δtᵢ)       # number of ice timesteps in the simulation
const subcyclenᵢ = Int(ΔtCpl/Δtᵢ)   # number of ice timesteps per coupling step
const outputΔt = 30minutes          # output frequency timestep
const outputnᵢ =  Int(outputΔt/Δtᵢ) # output frequency number for ice

# Physical constants and parameters:
const f = 1.4e-4                            # Coriolis parameter
const N² = (50f)^2                          # Stratification
const Ri = 1.0                              # Richardson number 
const Λ = sqrt(N²/Ri)                       # Background shear
const ρo = 1020.0                           # Ocean density
const hmean = 1.0                           # mean floe height
const Δh = 0.3                              # difference in floe heights - here all floes are the same height
const youngs_modulus = 1e5                  # Young's modulus for sea ice (a parameter that controls collision strength)
const z₀ᵢ = 6.0e-3                          # Ice roughness length (undeformed multiyear ice; e.g., McPhee, 2002; doi:10.1029/2000JC000633)
const κ = 0.41                              # Von Karman constant                 
const Cd_io = ( κ / log( (Δz/2)/z₀ᵢ) )^(2)  # Ice-ocean drag coefficient for the vertical grid resolution

# ## Ocean Model Setup
# The ocean model simulates of the Eady problem expanded around a background geostrophic shear with Ri= 1. 
# This approximately follows Listing 6 from [Wagner et al., 2025](https://arxiv.org/abs/2502.14148/v2).
#
# The ocean model setup follows similar steps as Subzero, with grid creation, model creation, simulation creation, 
# and output creation, along with steps for the background fields, and intial/boundary conditions. 
# For more information, see the [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) 
# documentation and examples.

# ### Grid creation
ocnGrid = RectilinearGrid( CPU(); size = (Nx, Ny, Nz),
                           x = (0,Lx), y = (0,Ly), z = (-Lz,0),
                           topology = (Periodic, Periodic, Bounded) )

# ### Background Fields
parameters = ( f=f, Λ=sqrt(N²/Ri) )
@inline U(x, y, z, t, p) = -p.Λ * z         # background velocity field
@inline B(x, y, z, t, p) = p.f * p.Λ * y    # background buoyancy field
background_fields = (u = BackgroundField(U; parameters),
                     b = BackgroundField(B; parameters) )
                     
# ### Boundary Conditions
# We initialize boundary conditions with empty fields that will be populated during coupling
Qu = Field{Center, Center, Nothing}(ocnGrid) 
Qv = Field{Center, Center, Nothing}(ocnGrid) 
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qu), ) 
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qv), )                     

# ### Model Creation
ocnModel = NonhydrostaticModel(; grid=ocnGrid, background_fields,
                               advection = WENO(), coriolis = FPlane(f=f),
                               tracers = :b, buoyancy = BuoyancyTracer(),
                               boundary_conditions = (; u=u_bcs, v=v_bcs ) )

# ### Initial Conditions 
ϵᵘ(x,y,z) = 1e-3 * randn() # initial random guassian noise seed
ϵᵇ(x,y,z) = 3e-4 * randn() # initial random guassian noise seed
bₒ(x,y,z) = N²*z + ϵᵇ(x,y,z) 
set!(ocnModel; b=bₒ, u=ϵᵘ, v=ϵᵘ )

# ### Simulation Creation
ocnSimulation = Oceananigans.Simulation(ocnModel, Δt=Δtₒ, stop_time=simTime )  

# ### Output Creation
dir = "coupled_ice_ocean"
mkpath(dir)
u,v,w = ocnModel.velocities
ζ = Field( (∂x(v)-∂y(u)) )
ocnSimulation.output_writers[:top] = JLD2Writer(  ocnModel, (; ζ); 
                                                  filename = dir * "/surface_vorticity",
                                                  schedule = TimeInterval(outputΔt),
                                                  overwrite_existing = true,
                                                  indices = (:, :, ocnGrid.Nz) )

# ### Add progress-update callback
function print_progress(sim)
    @printf("[%05.2f%%]   |   sim time: %8s   |   wall time: %8s   |   ocean CFL#: %3.2f\n",
            100 * (time(sim) / simTime),
            Oceananigans.prettytime( time(sim), false), 
            Oceananigans.prettytime( sim.run_wall_time, false ),
            AdvectiveCFL(sim.Δt)(sim.model)
            )
    flush(stdout)
    return nothing
end
add_callback!(ocnSimulation, print_progress, TimeInterval(0.025*simTime))



# ## Sea Ice Model Setup
# The setup for Subzero is the unchanged from normal except that in some cases we have to specify 
# the source for a function (to avoid conflicts with Oceananigans); for example, 
# `iceSimulation = Subzero.Simulation(...)` instead of `simulation = Simulation(...)`.


# ### Grid Creation
iceGrid = grid = RegRectilinearGrid(; x0 = 0.0, xf = Lx, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)

# ### Domain Creation
nboundary = PeriodicBoundary(North; grid=iceGrid)
sboundary = PeriodicBoundary(South; grid=iceGrid)
eboundary = PeriodicBoundary(East; grid=iceGrid)
wboundary = PeriodicBoundary(West; grid=iceGrid)
domain = Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary)

# ### Ocean and Atmosphere Instantiation
# We initialize the ocean field with zero-values that will be overwritten during coupling
atmos = Atmos(; grid=iceGrid, u = 0.0, v = 0.0, temp = 0.0 )
ocean = Ocean(; grid=iceGrid, u = 0.0, v = 0.0, temp = 0.0 )

# ### Floe Creation
floe_settings = FloeSettings( subfloe_point_generator = SubGridPointsGenerator(; Δg=Δgrid), # 1 internal point per gridpoint
                              maximum_ξ = 1 ) # increase maximum floe rotation rate
floe_generator = VoronoiTesselationFieldGenerator(; nfloes = 150, concentrations = [1/2], hmean, Δh)
floe_arr = initialize_floe_field(FT; generator = floe_generator, domain, rng = Xoshiro(1), floe_settings)


# ### Model Creation
iceModel = Model(; grid = iceGrid, ocean, atmos, domain, floes = floe_arr)

# ### Constants Creation
consts = Constants(; f = f, E = youngs_modulus, ρo = ρo, Cd_io = Cd_io, turnθ = 0 )

# ### Settings Creation
collision_settings = CollisionSettings( true, 1.0, 1.0 )
coupling_settings = CouplingSettings( true, subcyclenᵢ, 1, true )

# ### Output Creation
init_fn, floe_fn = "coupled_ice_ocean_init_state.jld2", "coupled_ice_ocean_floes.jld2"
initwriter = InitialStateOutputWriter(dir = dir, filename = init_fn, overwrite = true)
floewriter = FloeOutputWriter(outputnᵢ, dir = dir, filename = floe_fn, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# ### Simulation Creation
iceSimulation = Subzero.Simulation(; model = iceModel, consts, writers,  
                                    Δt = Int(Δtᵢ), nΔt = nΔtᵢ,
                                    floe_settings, collision_settings, coupling_settings,
                                    verbose = false )

# ## Coupling Setup
# We make use of the "Callback" functionality in Oceananigans to perform coupling between the two models.
# At each coupling timestep, surface fluxes and ocean currents are passed between the two models. Then, Subzero.jl
# is stepped forward in time for the number of ice-timesteps per coupling step.
# This is a serial coupling approach, where the two models are run one after another for each time segment.
#
# Functions to update surface fluxes and ocean currents

function updateSurfaceFluxes(sim) 
    Qu[0:Nx,0:Ny] .= -iceSimulation.model.ocean.τx./ρo  # m² s⁻²
    Qv[0:Nx,0:Ny] .= -iceSimulation.model.ocean.τy./ρo  # m² s⁻²
    return nothing
end

function updateOceanCurrents(sim) 
    Usurf = -Λ*(-Δz)/2 # background velocity at center of the first ocean grid cell below the surface
    iceSimulation.model.ocean.u .= ocnModel.velocities.u[0:Nx, 0:Ny, Nz] .+ Usurf
    iceSimulation.model.ocean.v .= ocnModel.velocities.v[0:Nx, 0:Ny, Nz] 
    return nothing
end

# Function to loop through the number of timesteps in the subcycle, running the model forward
tstepnᵢ = 0 # initial timestep number for sea ice
function runIceModel(sim)
    for _ in 1:subcyclenᵢ
        global tstepnᵢ
        Subzero.timestep_sim!( iceSimulation, tstepnᵢ )
        tstepnᵢ += 1
    end    
    return nothing    
end

# Adding Callbacks to Oceananigans Simulation
add_callback!(ocnSimulation, updateSurfaceFluxes, TimeInterval(ΔtCpl))
add_callback!(ocnSimulation, updateOceanCurrents, TimeInterval(ΔtCpl))
add_callback!(ocnSimulation, runIceModel, TimeInterval(ΔtCpl))


# ## Running the Coupled Simulation
# Because running Subzero occurs within an Oceananigans callback function, we run the full coupled model by calling
# the `run!` just for Oceananigans

Oceananigans.run!(ocnSimulation)

# ## Loading and Visualizing Results
# Instead of using the built-in `plot_sim` function, this example creates a custom visualization to 
# show both Subzero and Oceananigans outputs together.

# Load the results
floePolys = jldopen( joinpath(dir,"coupled_ice_ocean_floes.jld2"), "r")["poly"]
floeξ = jldopen( joinpath(dir,"coupled_ice_ocean_floes.jld2"), "r")["ξ"]
ζₛ_timeseries = FieldTimeSeries( joinpath(dir,"surface_vorticity.jld2"), "ζ" )
times = ζₛ_timeseries.times

# Use Makie's `Observable` to set up an animation of the data
floeKeys = keys(floePolys)
ni  = Observable( length(times) )
ζₛₙ = @lift ζₛ_timeseries[$ni]
floesₙ = @lift floePolys[ floeKeys[$ni] ]
floesζₙ = @lift 2.0.*floeξ[ floeKeys[$ni] ]
title = @lift @sprintf("t = %05.2f days", times[$ni]/day )

# Set up the figure and plot 2 panels showing ocean vorticity (both) and sea ice floes (left):
fig = Figure( )
ax = Axis( fig[1,1]; 
           limits = (0,Lx,0,Ly),
           xlabel = "x [km]", ylabel = "y [km]",
           xticks = ( (Lx/4).*(0:4), string.( (Lx/4kilometer).*(0:4) ) ),
           yticks = ( (Ly/4).*(0:4), string.( (Ly/4kilometer).*(0:4) ) ),
           )

heatmap!( ax, ζₛₙ; colormap = Reverse(:RdBu), colorrange = 2f.*(-1,1), )
poly!(ax, floesₙ; color=floesζₙ, strokecolor=:black, strokewidth=1, colormap = Reverse(:RdBu), colorrange = 2f.*(-1,1), )

ax = Axis( fig[1,2]; 
           limits = (0,Lx,0,Ly),
           xlabel = "x [km]", ylabel = "y [km]",
           xticks = ( (Lx/4).*(0:4), string.( (Lx/4kilometer).*(0:4) ) ),
           yticks = ( (Ly/4).*(0:4), string.( (Ly/4kilometer).*(0:4) ) ),
           )
hideydecorations!(ax,ticks=false)
hm = heatmap!( ax, ζₛₙ; colormap = Reverse(:RdBu), colorrange = 2f.*(-1,1), )
Colorbar( fig[1,3], hm; label = "ζ [s⁻¹]", ticks = ( f.*(-2:2), string.(-2:2).*"f" ) )

Label(fig[0, 1:2], title, fontsize=18, tellwidth=false)

rowsize!(fig.layout ,1, Aspect(1,1) )
display(fig)

# Record an animation:
@info "Making an animation..."
frames = 1:length(times)
record(fig, joinpath(dir,"oceananigans_coupled_example.mp4"), frames, framerate=24) do i
    ni[] = i
end
@info "Complete!"


# The results are shown in the animation below.  
# This animation shows the ocean surface vertical vorticity field, ζ, in both panels. 
# In the first panel (left), we overlay the ocean voriticy with the sea ice floes, each coloured by their own vorticity (ζ=2ξ). 
# In the second panel (right), we remove the ice floes, making it easier to see their imprints on the ocean.

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../coupled_ice_ocean/oceananigans_coupled_example.mp4" type="video/mp4">
# </video>
# ```