```@meta
EditURL = "../literate/tutorial.jl"
```

# Tutorial

This is the user tutorial for Subzero.jl. It will walk you through the typical workflow of
buidling a discrete-element model (DEM) with Subzero.jl as well as running and plotting your
simulatuion.

## Tutorial - copy-pasteable version

## Core ideas behind Subzero.jl simulations

The very first step of running a Subzero simulation is to bring the package into scope.

````@example tutorial
using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI
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
function shear_flow(Nx, Ny, min_u, max_u)
    increasing = true
    curr_u = min_u
    Δu = fld(Nx, 2)  # divide Nx by 2 and round down
    u_vals = zeros(Nx, Ny)
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
    return curr_u
end

u_vals = shear_flow(grid.Nx + 1, grid.Ny + 1, 0.0, 0.25)
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
atmos = Atmos(; grid, u = 5.0, v = 0.0, temp = 0.0)
````

Again since all of the fields are constant, we won't plot them, but you can, using the `heatmap`
function as shown above.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

