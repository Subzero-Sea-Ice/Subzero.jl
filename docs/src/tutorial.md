```@meta
EditURL = "../literate/tutorial.jl"
```

# Tutorial

This is the user tutorial for Subzero.jl. It will walk you through the typical workflow of
buidling a discrete-element model (DEM) with Subzero.jl as well as running and plotting your
simulatuion.

## Core ideas behind Subzero.jl simulations

The very first step of running a Subzero simulation is to bring the package into scope.

````@example tutorial
using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI
using Random, Logging
````

## Creating a Grid

Each Subzero model requires a grid object. The grid object defines the grid for the ocean
and atmosphere. Ocean and atmosphere vectors (like `u` and `v` velocities) and tracers
(like temperature) are recorded on these grid points and grid lines. The grid points are
then used for interpolation for coupling between the ice, ocean, and atmosphere.

All Subzero grid objects are concrete types of the abstract type [`AbstractRectilinearGrid`](@ref).
Right now, all grid objects must be Rectilinear. Currently the only implemented concrete type
is a [`RegRectilinearGrid`](@ref). If you are interested in implementing another type of
grid, see the [developer documentation]("devdocs.md").

Here, we will go ahead and create an instance of `RegRectilinearGrid`. We need to specify
the grid endpoints and either the number of grid cells in both directions, or the size of
the grid cells. Here we will specity the number of grid cells in the x-direction, `Nx`, and
in the y-direction, `Ny`.

````@example tutorial
grid = RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = 0.0, yf = 1e5, Nx = 20, Ny = 10)
````

We plot a dashed box around the grid so that we can see the that the grid matches the extent given.
We also place tick-marks at the desired grid cell lengths. Finally, set the plot's aspect ration to `2`
as the x-extent is two-times larger than the y-extent.

````@example tutorial
fig = Figure();
ax1 = Axis(fig[1, 1];  # set up axis tick marks to match grid cells
    title = "Grid Setup",
    xticks = range(grid.x0, grid.xf, 5), xminorticks = IntervalsBetween(5),
    xminorticksvisible = true, xminorgridvisible = true, xticklabelrotation = pi/4,
    yticks = range(grid.y0, grid.yf, 3), yminorticks = IntervalsBetween(5),
    yminorticksvisible = true, yminorgridvisible = true,
)
lines!(  # plot boundary of grid with a dashed line
    [grid.x0, grid.x0, grid.xf, grid.xf, grid.x0],  # point x-values
    [grid.y0, grid.yf, grid.yf, grid.y0, grid.y0];  # point y-values
    linestyle = :dash, linewidth = 3.0)
````

Resize grid to layout

````@example tutorial
colsize!(fig.layout, 1, Aspect(1, 2))
resize_to_layout!(fig)
fig  # display the figure
````

## Creating Boundaries

Next, each Subzero.jl model needs a `Domain`. A `Domain` defines the region of the grid that
the ice floes are allowed in, what happens to them when they reach the boundaries of that
region, and if there is any topography in the model, along with the ice, in that region.

Similarly to the `grid` above, the `Domain` will be rectilinear, defined by four boundaries,
one for each of the cardinal direction. You will be able to pass each of the cardinal
directions ([`North`](@ref), [`South`](@ref), [`East`](@ref), and [`West`](@ref)), defined as
types by Subzero, to the boundary constructors. Each boundary can have different behavior, allowing
the user to create a wide variety of domain behavior. Right now, four types of `AbstractBoundaries`
are implemented: [`OpenBoundary`](@ref), [`PeriodicBoundary`](@ref), [`CollisionBoundary`](@ref), and
[`MovingBoundary`](@ref).

In this example, we will use two `CollisionBoundary` walls and two `PeriodicBoundary` walls
to create a channel that the ice can infinitly flow through, from the east back to the west.
In the north and the south, the ice will collide with the boundaries, as if there was shoreline
preventing it from leaving the channel.

We will use the grid we made above to define the boundaries so that they exactly border the grid.

````@example tutorial
north_bound = CollisionBoundary(North; grid)
south_bound = CollisionBoundary(South; grid)
east_bound = PeriodicBoundary(East; grid)
west_bound = PeriodicBoundary(West; grid)
````

If we plot the polygons that are created within each boundary object, we can see that they border
the grid. These polygons are how we well when the ice floes are interacting with each of the
boundaries. We can also see that the boundaries overlap in the corners to ensure there is
a solid border around the grid. The `PeriodicBoundary` elements are in purple while the
`CollisionBoundary` elements are in teal.

````@example tutorial
poly!(   # plot each of the boundaries with a 50% transparent color so we can see the overlap
    [north_bound.poly, south_bound.poly, east_bound.poly, west_bound.poly];
    color = [(:purple, 0.5), (:purple, 0.5), (:teal, 0.5), (:teal, 0.5)],
)
ax1.title = "Grid and Boundary Setup"
fig  # display the figure
````

## Creating Topography
We then have the option to add in a [`TopographyField`](@ref Subzero.TopographyField), which is a collection of
[`TopographyElement`](@ref)s. If we want to add in topography field, we can create one using
the [`initialize_topography_field`](@ref) function. Here we will create two islands in the
channel. For simplcity, both will be triangles. I create the polygons that define the shape
of each island using [`GeoInterface`](https://github.com/JuliaGeo/GeoInterface.jl) and
defining the points with tuples:

````@example tutorial
island1 = GI.Polygon([[(-6e4, 7.5e4), (-4e4, 5e4), (-2.5e4, 7e4), (-6e4, 7.5e4)]])
island2 = GI.Polygon([[(5e4, 2.5e4), (5e4, 5.5e4), (7.5e4, 3e4), (5e4, 2.5e4)]])
topo_list = [island1, island2]
````

We can then pass these to `initialize_topography_field` with the `polys` keyword. We could
also have defined them just by their coordinates and passed in the coordiantes by the `coords`
keyword.

````@example tutorial
topo_field = initialize_topography_field(; polys = topo_list)
````

We can now plot the topography within the domain.

````@example tutorial
topo_color = RGBf(115/255, 93/255, 55/255)  # brown color for topography
poly!(topo_field.poly; color = topo_color) # plot the topography
ax1.title = "Grid and Domain Setup"
fig  # display the figure
````

## Creating a Domian
We now have all of the pieces needed to create a [`Domain`](@ref). We will combine the four
(4) boundaries we created, and the topography, into one `Domain` object. The collection of
boundaries and topography define where the floes can and cannot go in the simulation and add
in boundary behavior.We can do that as follows:

````@example tutorial
domain = Domain(; north = north_bound, south = south_bound, east = east_bound, west = west_bound, topography = topo_field)
````

!!! note
      We could have skipped adding the topography to the domain. That is an optional
      keyword and an empty `topography` field will be automatically created if the user does not
      provide one.

At this point, we have already plotted all of the `Domain` objects, so we will move in to adding
environmental forcing from the ocean and the atmosphere.

## Creating an Ocean
We can provide a ocean velocity and temperature fields to drive the floes within the simulatation.
By default, this is a static vector field and does not change throughout the simulation. However,
if you are interested in updating the code to make it easier to update the ocean/atmosphere fields
throughout the simulation please see issue [#111](@ref).

You can either provide a single constant value for each of the fields (and then also provide the `grid`
as input) or provide a user-defined field. If providing a user-defined field, the field should be of
size `(grid.Nx + 1, grid.Ny + 1)`. The rows represent velocities at the `x`-values of the grid
(from `x0` to `xf`), while the columns represent velocities at the the `y`-values of the grid
(from `y0` to `yf`). This makes it easy to index into the grid for values at point `(x, y)`.
If `x` is the `i`th x-value between `x0` to `xf` and `y` is the `j`the `j`th y-value between
`y0` to `yf` then you simply retrive the field values at index `[i, j]`. This does
mean that the fields provided don't "look" like the grid when output directly in the terminal
and would need to be transposed instead.

Here we will then define a `u`-velocity field that will impose a shear flow on the floes so
that ocean is flowing faster in the middle of the y-values (at 0.25 m/s) and then
symetrically slowing to 0m/s at the edges where the `CollisionBoundary` elements are. We will
then provide a constant `temp`erature of -1 C and a constant `v`-velocity of 0m/s.

````@example tutorial
function shear_flow_velocities(Nx, Ny, min_u, max_u)
    u_vals = fill(min_u, Nx, Ny)  # fill in field with min value
    increasing = true  # velocity values start by going up
    Δu = (max_u - min_u) / fld(Ny, 2)  # must reach max_u in the middle of the field
    curr_u = min_u + Δu  # first column already set to min_u
    for (i, col) in enumerate(eachcol(u_vals))
        (i == 1 || i == Ny) && continue  # edges already set to 0m/s
        col .= curr_u
        if increasing && curr_u == max_u  # reach max value and start decreasing velocity
            increasing = false
        end
        if increasing  # update velocity for next column
            curr_u += Δu
        else
            curr_u -= Δu
        end
    end
    return u_vals
end

u_vals = shear_flow_velocities(grid.Nx + 1, grid.Ny + 1, 0.0, 0.35)
ocean = Ocean(; u = u_vals, v = 0.0, temp = -1.0, grid)
````

We can now plot our ocean velocity values on an axis below our `grid` and `domain` setup.

````@example tutorial
ax2 = Axis(fig[2, 1]; title = "Ocean U-Velocities [m/s]", xticklabelrotation = pi/4)
xs = grid.x0:grid.Δx:grid.xf
ys = grid.y0:grid.Δy:grid.yf
u_hm = heatmap!(ax2, xs, ys, ocean.u)
Colorbar(fig[2, end+1], u_hm)
resize_to_layout!(fig)
fig
````

!!! note
      The [`heatmap`](https://docs.makie.org/v0.21/reference/plots/heatmap)
      already transposes the matrix input, assuming, as we do, that you index into the
      matrix with x-values for rows and y-values for columns. If you want to see your
      input as it would look on the grid, rather than manually transposing, just use the
      `heatmap` function.

Since the other ocean fields are constant values for each index, we won't plot these.

## Creating an Atmosphere
The atmosphere works almost iddntically to the `Ocean`. It also requires inputs for
`u`-velocities, `v`-velocities, and `temp`erature and these fields are considered static
throughtout the simulation by default. Again, if you are interested in changing this, please
see issue [#111](@ref).

Just like with the `Ocean` constructor, you can either provide the `Atmos` constructor a
single constant value for each of the fields (and then also provide the grid as input) or
provide a user-defined field.

Here, we will just provide constant values of 5m/s for the `u`-velocities, 0.0m/s for the
`v`-velocities and 0 C for the `temp`erature. If you are intersted in providing non-constant
fields, see the `Ocean` example above.

````@example tutorial
atmos = Atmos(; grid, u = 5.0, v = 0.0, temp = -2.0)
````

Again since all of the fields are constant, we won't plot them, but you can, using the `heatmap`
function as shown above.

## Creating the Floes

Within the simulation, there is a floe field, which is a list of [`Floe`](@ref) structs. For more information on these structs, see [`Floe`](@ref) documentation.
However, it is **not** recommended for users to create individual floes. Rather, users should use the `initialize_floe_field` function to create a field of floes
to start the simulation. Individual floes should only be created by developers when adding new simulation functionality.

### Creating a `AbstractFloeFieldGenerator`
To use [`initialize_floe_field`](@ref), you first need to create a concrete subtype of [`AbstractFloeFieldGenerator`](@ref). Right now, there are two implemented.
The first is the [`CoordinateListFieldGenerator`](@ref). This generator takes in a list of floe coordinates (along with other arguments) to create a list of floes.
The second is the [`VoronoiTesselationFieldGenerator`](@ref). This generator tesselates the domain with floes and then randomly keeps various floes to match the
user-provided floe concentration.

Here we will use the `VoronoiTesselationFieldGenerator` to generate a floe field. This generator object requires an input of the number of floes to be created (`nfloes`),
a vector/matrix of concentrations to split the field into (`concentrations`), a mean height for the generated floe field (`hmean`), and a maximum range of heights (`Δh`) such that
all floes have a height between `hmean - Δh` and `hmean + Δh`.

````@example tutorial
nfloes = 30
concentrations = [0.5 0.0]  # the left half of the domain is 85% packed and the right has no floes
hmean, Δh = 1.0, 0.25
generator = VoronoiTesselationFieldGenerator(; nfloes, concentrations, hmean, Δh)
````

### Creating `FloeSettings`

[`FloeSettings`](@ref) can be used to set things like the minimum allowed floe area, the maximum hieght, the density of ice, and the
`subfloe_point_generator`, which determines how to couple between the ice and ocean/atmosphere. There are also other fields within the
[`FloeSettings`](@ref)!! It is important to understand all of its fields.

!!! note
      For consistency, `FloeSettings` also must be passed into the [`Simulation`](@ref) struct.

Here we will just create a `FloeSettings` struct with just a few of the keyword options.

````@example tutorial
floe_settings = FloeSettings(;
  min_floe_area = 1e5,
  max_floe_height = 5,
  subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2),
)
````

### Calling `initialize_floe_field`

This generator is then the first keyword argument to the [`initialize_floe_field`](@ref) function. The other required keyword
is the simulation domain. After that, we have optional keywords, including but not limited to `supress_warnings`, `rng`, and `floe_settings`.
`supress_warnings` is a boolean that can turn off checks to ensure that all floe's have at least the minimum set area (see [`FloeSettings`](@ref) below)
and that all floe centroids are within the `domain`. `rng` is a random number generator, possibly seeded for reproducability, used to generate floe heights
and other random choices, like which floes to keep/remove to meet the requested concentrations.

````@example tutorial
floes = initialize_floe_field(; generator, domain, floe_settings, rng = Xoshiro(2))
````

We can add the floes to the domain above.

````@example tutorial
floe_color = RGBf(217/255, 226/255, 225/255)  # blue color for floes
poly!(ax1, floes.poly; color = floe_color) # plot the topography
ax1.title = "Grid, Domain, and Floe Setup"
fig  # display the figure
````

## Creating a Model

A `Model` combines all of the above components into the `Struct` that defines the physical aspects
of a Subzero simulation. There are a few constraints: the domain must fall within the grid, the ocean
and atmosphere must share the same grid (which can always be updated when the infrastructure for different
grids is in place), and all fields must share the same Float type (i.e. `Float64` or `Float32`). This Float
type will then be the Float type of the model, so it cannot be directly specified like it can when creating
the `Model` fields.

A model can be made as follows:

````@example tutorial
model = Model(; grid, domain, ocean, atmos, floes)
````

## Constants

We can then set a group of [`Constants`](@ref), which are a set of important physical parameters, like the ocean coriolis frequency and
the air density. All constants have a default value, so you only need to change the ones that specifically affect your simulation. All
constants and default values are listed in the [`Constants`](@ref) documentation.

Here, I will just use the default values:

````@example tutorial
consts = Constants()
````

## Settings

In addition to the [`FloeSettings`](@ref) discussed above, there are also settings for all of the other physical processes that can happen
during a Subzero run. Here, I will just use all default settings except for the `FractureSettings`, but each of the following settings have documentation
detailing how to turn features on/off and tune their behavior.

The list of existing settings is:
- [`FloeSettings`](@ref)
- [`CouplingSettings`](@ref)
- [`CollisionSettings`](@ref)
- [`FractureSettings`](@ref)
- [`SimplificationSettings`](@ref)
- [`RidgeRaftSettings`](@ref)
- [`WeldSettings`](@ref)

````@example tutorial
fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floes),
        Δt = 500,
        npieces = 3,
        deform_on = false,
)
````

## Output Writers

You also need to create a set of output writers to record the data from your simulation. You can add four types of output writers, and as many of each
type as you would like. The four types are as follows: [`InitialStateOutputWriter`](@ref), [`CheckpointOutputWriter`](@ref),[ `FloeOutputWriter`](@ref), and
[`GridOutputWriter`](@ref). When any of these objects are created, the file that they will write to is also created automatically. You can read more about each
of these types of output writers within the API guide. You must then combine all of your output writers into an [`OutputWriters`](@ref) struct. You can simply
pass each output writer created as an argument to the `OutputWriters` constructor.

Here, we will create an [`InitialStateOutputWriter`](@ref), a [`CheckpointOutputWriter`](@ref) (that outputs every 1000 timesteos), and a
[ `FloeOutputWriter`](@ref) (that outputs every 50 timesteps) and then combine them into a [`OutputWriter`](@ref).

````@example tutorial
dir = "tutorial"
init_fn, checkpoint_fn, floe_fn = "tutorial_init_state.jld2", "tutorial_checkpoint.jld2", "tutorial_floes.jld2"
````

We first make the [`InitialStateOutputWriter`](@ref):

````@example tutorial
initwriter = InitialStateOutputWriter(; dir = dir, filename = init_fn, overwrite = true)
````

We then make the [`CheckpointOutputWriter`](@ref):

````@example tutorial
checkpointer = CheckpointOutputWriter(1000; dir = dir, filename = checkpoint_fn, overwrite = true)
````

Finally, we make a [`FloeOutputWriter`](@ref):

````@example tutorial
floewriter = FloeOutputWriter(50; dir = dir, filename = floe_fn, overwrite = true)
````

We can then combine these into an `OutputWriters` object:

````@example tutorial
writers = OutputWriters(initwriter, checkpointer, floewriter)
````

## Create a Simulation

At this point, you are ready to make a `Simulation`. A `Simulation` has quite a few fields, all of which are talked about above in more detail.
Since there are so many fields, they are keyword defined, so you must provide a keyword when creating the struct. The only necessary arguments
are the `model` (to specify the simulation physical setup) and `floe_settings` (to ensure they match the argument used when creatting the `floes`).
See the [`Simulation`](@ref) documentation for a full list of keyword arguments and their default values.

In this tutorial, we create a simulation using all of the structs we created above. We will also need to set the simulation timestep `Δt` and the
total number of timesteps to run for `nΔt`.

````@example tutorial
sim = Simulation(;
    model = model,
    consts = consts,
    Δt = 10, # timestep of 5 seconds
    nΔt = 20000, # run for 20,000 timesteps
    floe_settings = floe_settings,
    fracture_settings = fracture_settings,
    writers = writers,
)
````

!!! note
      Subzero also has a custom logger, [`SubzeroLogger`](@ref), which logs all of the info and warning messages that Subzero throws over
      the course of a simularion into a log file. It will only output each unique message `messages_per_tstep` times, which can be passed
      to the `run!` function. You can also create your own `SubzeroLogger`. However, here, for simplicity of the tutorial, I will turn off
      all logging at or below `Info`.

````@example tutorial
Logging.disable_logging(Logging.Info)
````

## Running the Simulation
You can now use the [`run!`](@ref) function to run the simulation.

````@example tutorial
run!(sim)
````

!!! note
      If you wish to couple to Oceananigans, you will need to run each model timestep by timestep
      and pass the needed fields back and forth. You can run a single timestep of the simulation
      using the [`timestep_sim!`](@ref) function.

If you run your simulation in multiple parts and need to re-start your simulation from files,
the [`restart!`](@ref) function will be a good place to start. However, note that it is quite simple
and users may need to write their own restart function if they want more complex behavior.

## Plotting the Simulation

If your simulation has both a [`FloeOutputWriter`](@ref) and an [`InitialStateOutputWriter`](@ref),
you can use the built in plotting function to make an MP4 file with each frame as a timestep saved
by the `FloeOutputWriter`. This plotting function `plot_sim` is quite simple and just meant to get you started.
You may need to add more complex plotting code to suit your needs.

````@example tutorial
plot_sim(joinpath(dir, floe_fn), joinpath(dir, init_fn), sim.Δt, joinpath(dir, "tutorial.mp4"))
````

```@raw html
<video width="auto" controls autoplay loop>
<source src="../tutorial/tutorial.mp4" type="video/mp4">
</video>
```

## Reproducibility

The simulations are currently completly reproducible when run single-threaded. The simulation takes a `rng` argument, and if it is provided with a seeded random number generator,
then two simulations with the same set of starting floes will produce the exact same results. Note that the same set of floes can be reproduced by providing a seeded random number
generator to the floe creation functions. Simulation reproducability has not been achieved for multi-threaded runs.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

