# How to Create a Simulation: 

You have a lot of flexibility in designing a simulation within Subzero. First, you need to create all the pieces of your model, including the grid, the ocean, the atmosphere, the domain, and the floes. After than, you can use your model to create a simulation, where you will also specify physical constants and other runtime parameters.

## Contents
* [Building a Model](#building-a-model)
    * [Model Calculations Type](#model-calculations-type)
    * [Grid](#grid)
    * [Ocean](#ocean)
    * [Atmosphere](#atmosphere)
    * [Domain](#domain)
    * [Floes](#floes)
    * [Making the Model](#making-the-model)
* [Building a Simulation](#building-a-simulation)
    * [Constants](#constants)
    * [Physical Process Settings](#physical-process-settings)
    * [Timesteps](#timesteps)
    * [Output Writers](#output-writers)
    * [Reproducibility](#reproducibility)
    * [Creating the Simulation](#creating-the-simulation)

## Building a Model
### Model Calculations Type: 

Simulations are designed to run with either Float64 or Float32 calculations. However, as of now, only Float64 is tested and supported. When you are creating a model and a simulation, you need to create all elements in the same type that you are trying to run the model in. To acomplish this, there is an optional first argument specifying either Float64 or Float32 for all functions that create model/simulation objects. However, Float64 is the default and the type argument can be dropped and all elements will be made with Float64.

To explcitly set this optional argumentAt the beginning of each example file, you will see the line: 

```julia 
const FT = Float64 
``` 
It is then used in each of the constructors, as can be seen below. Note again that this argument can be dropped and it will default to Float64.

### Grid: 
The first thing that you need to create in your grid. For now, the ocean and atmosphere must be on the same grid and the only grid is a RegularRectilinearGrid, where every grid cell is a rectangle of the same size. You can create a regular rectilinear grid one of two ways. Both methods require specifying the minimum and maximum x and y points in meters. 

The first method has you provide the size of each grid cell with ∆x and ∆y arguments. That can be seen here: 

```julia 
grid = RegRectilinearGrid( 
  FT, 
  (-1e5, 1e5), # x bounds
  (0.0, 1e5),  # y bounds
  2e4,         # ∆x
  1e4,         # ∆y
) 
``` 

Which will create a grid that has a width of 2e5m from –1e5m to 1e5m and each grid cell was a width of 2e4 and a height of 1e4. This would then create a grid with 100 cells, 10 cells across and 10 cells wide. Note that if you provide grid cells dimensions that don’t divide the grid width and height evenly, the grid will be truncated to the nearest grid cell. 

The other way to create this same grid is to use: 

```julia 
grid = RegRectilinearGrid(
   FT,
   10,  # Nx
   10,  # Ny
   (-1e5, 1e5),  # x bounds
   (0.0, 1e5),   # y bounds
) 
``` 
by specifying the desired grid dimensions. Here you are guaranteed to always get a 10 by 10 grid, but do not have as fine-grained control over the size of your cells.  

A grid object has the following fields:
Nx, Ny, x0, xf, y0, yf, Δx, and Δy
Ny is the number of rows and Ny is the number of columns within the grid as defined by each grid cell. x0 and xf are the minimum and maximum x-bounds, and y0 and yf are the minimum and maximum y-bounds. Finally, Δx and Δy are the number of grid cells in the x and y directions.

Additionally, RegRectilinearGrid is a concrete subtype of AbstractRectilinearGrid. More concrete subtypes may be added in the future.  

### Ocean: 
The ocean here represents 2D vector fields of the surface layer of the ocean. It includes the following: u-velocity, v-velocity, temperature, and stresses in the x and y direction, the amount of area per cell covered in ice, and the `hflx_factor` per cell, which allows local heat flux calculations based on the difference in atmosphere and ocean temperature in each cell. If you want to run Subzero without coupling, the ocean will be a set of prescribed fields. If you couple Subzero to Oceananigans, Oceananigans will provide the velocity and temperature fields, and Subzero will provide oceananigans with stress fields from the ice and atmosphere on the top layer of the ocean.

There are several ways to initialize an ocean using Subzero. If you want to create uniform velocity and temperature fields, those constant values can simply be provided. This is useful for testing as a quick way to create a simple environment. This can be done as follows:

```
ocean = Ocean(FT, grid, -0.3, 0.2, 0.0)
```
to create an ocean with a -0.3m/s zonal velocity, a 0.2m/s meridional velocity, and a 0.0°C temperature field. The other ocean fields, such as the ocean stresses, will be calculated through out the simulation. 

To create more interesting fields, the ocean fields must be defined and explicitly provided to the Ocean constructor. To create an equivalent ocean as above, the following syntax can be used:
```
ocean = Ocean(
   FT,
   fill(-0.3, 11, 11), # ocean values are stored on grid lines
   fill(0.2, 11, 11),
   zeros(FT, 11, 11)
)

Note that the ocean velocity matricies are x values by y value (i.e. the x values are rows and the y values are columns) for ease of indexing and to match with Oceananigans. 
```

### Atmosphere: 
The atmosphere is very similar to the ocean in that is also represents a collection of 2D vector fields. We have not yet coupled Subzero with any atmosphere model. Therefore, for now, the atmosphere will be perscribed. It also has a u-velocity, v-velocity, and temperature field and can be created identically to the ocean. The two methods for creating the atmosphere are shown below:
```
atmos = Atmos(FT, grid, -0.3, 0.2, 0.0)

atmos = Atmos(
   FT,
   fill(-0.3, 11, 11), # atmosphere values are stored on grid lines
   fill(0.2, 11, 11),
   zeros(FT, 11, 11)
)
```

### Domain: 
The domain includes both the boundaries of the simulation and the topographic features. It defines the areas where the floes cannot go. Before you can define a domain, you need to define each boundary wall. As of now, only rectangular boundaries around the edge of the grid are allowed so we need a north, east, south, and west boundary wall. We have four types of boundary wall: open, collision, periodic, and moving.

With an open wall, if a floe overlaps with the boundary wall at all, the floe is removed from the simulation. With a collision wall, floes collide with the wall, and it is calculated similarly to a collision with another floe, except that since the wall cannot break, or be moved, it has an idealized force factor. With a periodic wall, if a floe overlaps with the wall, a copy of that floe is created passing back into the domain through the opposite wall. These two floes are equivalent and are linked. We call the copy a “ghost floe.” So, for example, if a floe is pushed partially out of the domain through the south wall by a current, a copy of the floe will be created, re-entering through the north wall. If one wall in the domain is periodic, its opposite wall must also be periodic, that is north-south and east-west. Moving boundaries are walls that can move with a constant velocity in either the x or y-direction, causing either a shear or compressive stress. For example, if the west wall is a moving wall with a x-velocity of 0.1m/s it will move with a 0.1m/s u velocity and a 0m/s v velocity towards the center of the domain. These types of walls are for increasing pressure on the ice to investigate stress and strain on the floes.  

Here is an example of creating a set of boundary walls using the grid from above: 

```julia 
nboundary = PeriodicBoundary(North, grid) 
sboundary = PeriodicBoundary(South, grid) 
eboundary = CollisionBoundary(East, grid) 
wboundary = OpenBoundary(West, grid) 
``` 

Once we have defined our walls, we can create a domain as follows: 

```julia 
domain = Domain(nboundary, sboundary, eboundary, wboundary)
``` 

Here is an example of how to create a moving wall as well that moves towards the center of the domain with u velocity of 0m/s and a v velocity of -0.1m/s:
```julia 
move_boundary = MovingBoundary(North, grid, 0.0, -0.1)
```


However, if we want topography, we can also add it to the domain. To do that, we would create one, or more, TopographyElements. Let us consider two simple square islands: 

```julia 
island1 = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
island2 = [[[8e4, 4e4], [8e4, 4.5e4], [8.5e4, 4.5e4], [8.5e4, 4e4], [8e4, 4e4]]] 
topo_arr = initialize_topography_field(FT, [island1, island2]) 
``` 

We can then add an additional argument to the domain when it is created:  

```julia 
domain = Domain(nboundary, sboundary, eboundary, wboundary, topo_arr) 
``` 

These are the two ways to define a domain in Subzero.  

For now, you will need to define ocean velocities that go around your topography, as Subzero does not change the given ocean velocities. This holds true on coupled scenarios as well. You need to make sure that if you do couple Subzero to Oceananigans that you set the ocean height to 0m at these locations so that the current flows around the topography. 

### Floes

Floes are quite complex objects as they need a lot fields. Here we will talk about a floe struct's fields, as well as how to create a configuration of floes to start your simulation.

#### Floe Struct Fields
A floe's fields can be broken down into several catagories. We will go through each catagory and describe the fields within in briefly in table-form.

The first catagory is **physical properties**. These have to do with the floe's physical shape. 

Before listing the fields, one important thing to know is that a floe's coordinates are represented by a `PolyVec`, which is a shorthand for a vector of a vector of a vector of floats. This sounds complicated, but it is simply a way of representing a polygon's coordinates. A Polygon's coordinates are of the form below, where the xy-coordinates are the exterior border of the floe and the wz-coordinates, or any other following sets of coordinates, describe holes within the floe:

```julia
coords = [
  [[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]],  # Exterior vertices of the polygon represented as a list of cartesian points
  [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]],  # Interior holes of the polygon represented as a list of cartesian points
  ...,  # Additional holes within the polygon represented as a list of cartesian points
 ]
 ```
 We will use the term `PolyVec` to describe this form of coordiantes and you will see it in the code if you take a look at the code base. It is also the form that floe coordinates are saved in output files.
 
| Physical Fields| Meaning                            | Type          |
| -------------- | ---------------------------------- | ------------- |
| centroid       | floe's centroid                    | Float64 or Float32|
| coords         | floe's coordinates                 | PolyVec of Float64 or Float32 |
| height         | floe's height in [m]                 | Float64 or Float32|
| area           | floe's area in [m^2]                 | Float64 or Float32|
| mass           | floe's mass in [kg]                  | Float64 or Float32|
| rmax           | floe's maximum radius, the maximum <br> distance from centroid to vertex in [m] | Float64 or Float32|
| moment         | floe's mass moment of intertia in [kg m^2]    | Float64 or Float32|
| angles         | list of floe's vertex angles in [degrees] | Vector of Float64 or Float32|

The second catagory is **sub-floe points**. These are used for interpolation of the ocean and atmosphere onto the floe. They are a list of points within the floe. There are two ways they can be generated. One is for them are randomly generated, with the user providing an initial target number of points. These are generated with a monte carlo generator (explained below). The other way is for the points to be on a sub-grid within the floe. There are benefits and drawbacks to both strategies. 

| Sub-floe Point Fields| Meaning                                 | Type                        |
| ----------------- | --------------------------------------- | --------------------------- |
| x_subfloe_points   | floe's sub-floe points x-coordinates | Vector of Float64 or Float32|
| y_subfloe_points   | floe's sub-floe points points y-coordinates | Vector of Float64 or Float32|

The third catagory is **velocities and orientations**. Floe's have both linear and angular velocity and keep track of the angle that they have rotated since the begining of the simulation.
| Movement Fields| Meaning                         | Type               |
| -------------- | ------------------------------- | ------------------ |
| u              | floe's x-velocity in [m/s]        | Float64 or Float32 |
| v              | floe's x-velocity in [m/s]        | Float64 or Float32 |
| ξ              | floe's angular velocity in [rad/s]| Float64 or Float32 |
| α              | rotation from starting position<br> in [rad]| Float64 or Float32 |

The fourth catagory is **status**. These fields hold logistical information about each floe and through which process it originated.
| Status Fields  | Meaning                         | Type               |
| -------------- | ------------------------------- | ------------------ |
| status         | if the floe is still active in the simulation        | Subzero.Status (see below)|
| id             | unique floe id for tracking the floe throughout the simulation | Int |
| ghost_id       | if floe is not a ghost, `ghost_id = 0`, else it is in `[1, 4]`<br> as each floe can have up to 4 ghosts| Int |
| parent_id    | if floe is created from a fracture or the fusion of two floes, `parent_id` is a list of <br> the original floes' `id`, else it is emtpy | Vector of Ints |
| ghosts         | indices of floe's ghost floes within the floe list| Vector of Ints |

A Status object has two fields: a tag and a fuse index list. There are currently three different tags: `active`, `remove`, and `fuse`. If a floe is `active`, it will continue in the simulation at the end of a timestep. If a floe's tag is `remove`, it will be removed at the end of the timestep. This ususally happens if a floe exits the domain, or becomes unstable for some reason. If a floe is marked at `fuse`, this means that is is overlapping with another floe by more than the user defined maximum overlap percent (see [Physical Process Settings](#physical-process-settings) for more information on this maximum fraction value. If a floe is marked for fusion, the index of the floe it is supposed to fuse with will be listed in the `fuse_idx` list.  

The fifth catagory is **forces and collisions**. These fields hold information about the forces on each floe and the collisions it has been in.
| Force Fields      | Meaning                                    | Type                        |
| ----------------- | ------------------------------------------ | --------------------------- |
| fxOA              | x-force on floe from ocean and atmosphere in [N] | Float64 or Float32|
| fyOA              | y-force on floe from ocean and atmosphere in [N] | Float64 or Float32|
| trqOA             | torque on floe from ocean and atmosphere in [N m]| Float64 or Float32|
| hflx_factor       | coefficent of floe height to get heat flux directly <br> under floe in [W/m^3]| Float64 or Float32|
| overarea          | total overlap of floe from collisions in [m^2]   | Float64 or Float32|
| collision_force   | forces on floe from collisions in [N]            | Float64 or Float32|
| collision_trq     | torque on floe from collisions in [N m]          | Float64 or Float32|
| interactions      | each row holds one collision's information, see below for more information | `n`x7 Matrix of Float64 or Float32 <br> where `n` is the number of collisions|
| stress_accum      | stress accumulated over the floe over past timesteps given StressCalculator, where it is of the form [xx yx; xy yy] | 2x2 Matrix{AbstractFloat}|
| stress_instant    | instantaneous stress on floe in current timestep | 2x2 Matrix{AbstractFloat} 
| strain            | strain on floe where it is of the form [ux vx; uy vy] | 2x2 Matrix of Float64 or Float32|

The `interactions` field is a matrix where every row is a different collision that the floe has experienced in the given timestep. There are then seven columns, which are as follows:
- `floeidx`, which is the index of the floe that the current floe collided with
- `xforce`, which is the force in the x-direction caused by the collision
- `yforce`, which is the force in the y-direction caused by the collision
- `xpoint`, which is the x-coordinate of the collision point, which is the x-centroid of the overlap between the floes
- `ypoint`, which is the y-coordinate of the collision point, which is the y-centroid of the overlap between the floes
- `torque`, which is the torque caused by the collision
- `overlap`, which is the overlap area between the two floes in the collision. 
You can use these column names to access columns of `interactions`. For example: `floe.interactions[:, xforce]` gives the list of x-force values for all collisions a given floe was involved in as that timestep.

The fifth catagory is **previous values**.
| Previous Value Fields | Meaning                                       | Type               |
| --------------------- | --------------------------------------------- | -------------------|
| p_dxdt                | previous timestep x-velocity (u) in [m/s]         | Float64 or Float32|
| p_dydt                | previous timestep y-velocity (v) in [m/s]         | Float64 or Float32|
| p_dudt                | previous timestep x-acceleration in [m/s^2]   | Float64 or Float32|
| p_dvdt                | previous timestep y-acceleration in [m/s^2]   | Float64 or Float32|
| p_dαdt                | previous timestep angular-velocity in [rad/s] | Float64 or Float32|
| p_dξdt                | previous timestep time angular acceleration in [rad/s^2] | Float64 or Float32|

#### Floe Settings

When you create a floe or a set of floes, you have the option to create a floe settings object. This set of settings controls certian floe fields and calculations.

The following fields are part of the floe settings (with default values):
  - ρi: floe's density (920.0 g/L)
  - min_floe_area: minimum floe area (1e6 m^2)
  - min_floe_height: minimum floe height (0.1 m)
  - max_floe_height: maximum floe height (10.0 m)
  - min_aspect_ratio: minimum ratio between floe x-length and y-length by maximum coordiante values (0.05)
  - maximum_ξ: the absolute maximum rotational velocity a floe can reach before it is capped at maximum_ξ (1e-5 rad/s)
  - subfloe_point_generator: generates floe's subfloe points (`MonteCarloPointsGenerator()`)
  - stress_calculator: generates the calculator for stress ('DecayAreaScaledCalculator()')

If any of the minimum values are exceeded, a floe is removed in the course of the simulation. If any of the maximum values are reached, the value is capped at the given value.

There are two types of subfloe point generators.

The first is a `MonteCarloPointsGenerator`. This generates random points within the floe. You can create a MonteCarloPoint generator with three fields: `npoints` (default 1000), `ntries` (default 100), and `err` (default 0.1). `npoints` is the number of points to attempt to generate, `ntries` is the number of tries to generate a set of points that meets the acceptable error, and `err` is the percent of floe are that can not be covered by monte carlo points for it to be a valid set of subfloe points. The user will not end up with `npoints` monte carlo points. These are the number of points generated in a bounding box around the floe. However, every point that is outside of the floe will be removed. The monte carlo points are repeatedly generated until a set is created with less that `err`. If a set cannot be determined in `ntries` tries, the floe will be marked for removal using the `status` field (see below).

The second type of subfloe point generator is the `SubGridPointsGenerator`. This generator places points along a grid within the floe. The user can define how fine that grid should be in comparison with the model's grid. A `SubGridPointsGenerator` takes in two arguemnts: the model's `grid` and `npoint_per_cell`, which defines how many subfloe points the user wants within the model's grid cell in botht the x and y direction (i.e. `npoint_per_cell = 3` will give 9 points in a grid cell, three in both x and y in a grid pattern).

Both of these have different benefits. `MonteCarloPointsGenerator` guarentee that all floes have a somewhat similar number of subfloe points. However, since these points are randomly placed, they are not neccesarily evenly spread out. There may not even be one point per model grid cell, which causes errors when two-way coupling as stress from ice to ocean is calcualted with subfloe points. 

On the other hand with a `SubGridPointsGenerator`, each floe has a number of points proportional to its area. However, these points are mainly evenly spaced and it is guarenteed that there is at least one per every model grid cell. Therefore, you must use `SubGridPointsGenerator` when two-way coupling. 

You can make a floe settings object as follows:
```julia
floe_settings = FloeSettings(
  min_floe_area = 1e5,
  max_floe_height = 5,
  subfloe_point_generator = SubGridPointsGenerator(grid, 2),
  stress_calculator = DecayAreaScaledCalculator(50),
)
 ```
 Any fields that aren't specified are assigned their default value.

#### Construct Individual Floes

You can create one floe at a time using floe constructors that will set initial values for all of these fields depending on your inputs.

Here is an example of using the PolyVec coordinates constructor, assume we have already created a PolyVec called `coords`:
```julia
coords = [[[3e4, 1e4], [3e4, 1.5e4], [3.5e4, 1.5e4], [3.5e4, 1e4], [3e4, 1e4]]]
hmean = 0.25
Δh = 0.1
floe = Floe(
    FT,
    coords,
    hmean,  # Floe height will be between 0.15 - 0.35
    Δh;  # Δh is the maximum difference between hmean and actual floe height
    floe_settings = floe_settings,
    u = 0.0,
    v = 0.0,
    ξ = 0.0,
    rng = Xoshiro(1), # seed of 1
)
```
Only `coords`, `hmean`, and `Δh` are neccesary arguments. The rest aer optional, and the default values are the values shown in the code snippit above.

However, it is not recomended that you manually create each floe. It is recomended that you use the `initialize_floe_field` functions instead to create your simulation's starting configuration of floes.
 
#### Initial Floe Configuration
It is recomeneded that you use the `initialize_floe_field` to create your starting configuration on floes. There are two ways to use this function. The first way is to provide a list of `PolyVecs` representing a list of the coordinates of all of the floes you want in the initial state. This will initialize all of the given polygons specified as floes. The other way is to provide a number of floes and a concentration over a specific area. This will create a starting floe field using voronoi tesselation that aims to achieve the requested number of floes and concentrations.

Both of these functions share most arguments. They are as follows (with default values if they exist):
- `domain`, which is the model's domain so that the floes can be fit into the open space using polygon intersections/differences
- `hmean`, which is the mean height of all floes created
- `Δh`, which is the maximum potential height difference from `hmean` between floes
- `floe_settings = FloeSettings()`, which specifies many values needed to create floes, as detailed above
- `rng = Xoshiro()`, which is a random number generator so that floe creation is reproducible if a seeded random number generator is provided

Note that all arguments with default values are optional keyword arguments.

Here is an example of creating a small floe field using the version of `initialize_floe_field` that takes in lists of `PolyVec`s. For a real simulation, you would probably generate a list of coordinates and read them in from a file but we create these by hand here for simplicity. Assume we have already created a `domain`.
```julia
floe1 = [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
floe2 = [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]
floe_field = initialize_floe_field(
  FT,
  [floe1, floe2],
  domain,
  0.25,  # mean height of 0.25
  0.0;  # all floes will be the same height
  rng = Xoshiro(1),
  floe_settings = floe_settings,
)
```

Now here is an an example of creating a large floe field using the version of `initialize_floe_field` that uses Voronoi tesselation. Again assume we have already created a `domain`.
```julia
floe_arr = initialize_floe_field(
    FT,
    100,  # attempt to initialize 100 floes
    [1.0; 0.0],  # the top half of the domain is fully packed and the bottom has no floes
    domain,
    0.25,  # mean height of 0.25
    0.10;  # floe heights will range from 0.15-0.35
    floe_settings = floe_settings,
    rng = Xoshiro(1),
)
```
We now focus on the first two arguments. The first is the number of floes to attempt to create with Voronoi tesselation. We are not guarenteed to get exactly that number. It depends on the amount of open space in the domain and the generation of random seed points. For example, if the domain is filled with lots of topography and islands, it will be more difficult to hit the exact number of floes requested. However, it will be in the ballpark. The second argument is the concentrations, which is a matrix. We can split the domain into quadrents that are the same shape at matrix and then request concentrations of ice in each of those quadrents equal to the corresponding value in the concentrations matrix. The other arguments are the same as in the floe coordinate version on the function.

### Making the Model
Once you have made all of the above components, you are now able to make a model. You will do that as follows:
```julia
model = Model(grid, ocean, atmos, domain, floe_arr)
```
This model object now holds all of the previously created components. You will get an error if you did not make all of your componenets with the same `FT` type. If you did make all objects with the same value of `FT`, either `Float64` or `Float32`, then the model now has that type.

## Building a Simulation
Your simulation will hold your model, as well as some runtime and physical parameters.

### Constants
Running the simulation requires known physical parameters, such as drag coefficents, densities, and the coriolis parameter. If you want to change any of the default values, you will need to create a `Constants` object. You can then keyword define any of the values that you want to change. All others will stay with their default value. For example, below we are creating a constants object that has a coefficent of friction set to 0, to run the simulation without friction. 

```julia
consts = Constants(μ = 0.0)
```
The list of constants and their default values are shown below:
| Constants     |    Meaning                         | Default Value |
| ------------- | ---------------------------------- | ------------- |
| ρo            | ocean densisty                     | 1027.0 kg/m^3 |
| ρa            | air densisty                       | 1.2 kg/m^3    |
| Cd_io         | ice-ocean drag coefficent          | 3e-3          |
| Cd_ia         | ice-air drag coefficent            | 1e-3          |
| Cd_ao         | air-ocean drag coefficent          | 1.25e-3       |
| f             | ocean coriolis frequency           | 1.4e-4 rad/s  |
| turnθ         | ocean turning angle                | 15π/180 rad   |
| L             | latent heat of freezing            | 2.93e5 J/kg   |
| k             | thermal conductivity of surface ice| 2.14 W/(mK)   |
| ν             | Poisson's ratio                    | 0.3           |
| μ             | Coefficent of friction             | 0.2           |
| E             | Young's Modulus                    | 6e6 N/m2      |

In particular, Young's Modulus is usually calculated using the total floe area after floe initialization in the original Subzero code:
```julia
E = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
```

### Physical Process Settings
Subzero allows you to turn on and off various physical processes, as well as change the settings for these physical processes. The physical processes availible are: coupling, collisions, fractures, and simplfication. We will be adding corner fracturing, packing, rafting, ridging, and welding. For each of these physical processes, you can create a settings object, and change the parameters within that object to change the physical process to meet your simulation needs. 
#### Coupling Settings
`CouplingSettings` changes how Subzero uses the ocean and atmosphere two-dimensional vector fields within the `Ocean` and `Atmos` structs, whether those values are static user-provided values or values that are updated every timesteps through coupling with Oceananigans. The availible settings are:
- coupling_on, which turns this process on and off
- ∆t, which sets the number of timesteps between this process running
- ∆d, which sets the number of ocean/atmosphere grid cells around a floe to consider when interpolating ocean/atmosphere forcings
- two_way_coupling_on, which turns on and off the calculation of the ice and atmosphere effects on the ocean and stores output in the ocean's stress fields.

Here is an example of creating your own coupling settings, using the default values:
```julia
couple_settings = CouplingSettings(
  coupling_on = true,
  Δt = 10,
  Δd = 1,
  two_way_coupling_on = false,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Collision Settings
`CollisionSettings` changes Subzero's floe-floe and floe-domain collisions. The availible settings are:
- collisions_on, which turns this process on and off
- floe_floe_max_overlap, which sets the maximum fraction of floe-floe overlap before fusing the two floes together 
- floe_domain_max_overlap, which sets the maximum fraction for floe-domain overlap before removing the floe from the simulation

Here is an example of creating your own collision settings, using the default values:
```julia
collision_settings = CollisionSettings(
  FT,
  collisions_on = true,
  floe_floe_max_overlap = 0.55,
  floe_domain_max_overlap = 0.75,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.
#### Fracture Settings
`FractureSettings` changes Subzero's floe fractures. The availible settings are:
- fractures_on, which turns this process on and off
- criteria, which sets the rules for which floes fractures (see more below)
- Δt, which sets the number of timesteps between this process running
- deform_on, which turns on a sub-process where a floe is deformed around the floe that is overlaps with the most during that timestep before it is fractured
- npieces, which sets the number of pieces a floe should try to fracture into (note that this number might not be met depending on floe shape)

Here is an example of creating your own fracture settings, using the default values:
```julia
fracture_settings = FractureSettings(
  fractures_on = false,
  criteria = NoFracture(),
  Δt = 0,
  deform_on = false,
  npieces = 3,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

If you do want to turn fractures on, you will need to create a criteria that is not the `NoFracture()` criteria. Right now, there is are two other fracture criteria.

One is Hibler's elliptical yield curve. To learn more about this criteria, see the 1979 paper "A Dynamic Thermodynamic Sea Ice Model." This criteria uses the mean floe height and take in two parameters, `pstar` and `c`, which can be tuned to get the fracture behavior that you want. The `HilberYieldCurve` object then has both of these parameters as fields as well as the `vertices` of the yield curve. A simple way to create this criteria is as follows:
```julia
pstar = 2.25e5
c = 20.0
criteria = HiblerYieldCurve(FT, floe_arr, pstar, c)
```

The other is Mohr's Cone yield curve. To learn more about this critera, see the 2009 paper "Coulombic faulting from the grain scale to the geophysical scale: lessons from ice." This criteria takes in three parameters that define the shape of the cone. They are as follows:
- q: based on the coefficient of internal friction (µi) by (μi^2 + 1)^(1/2) + μi^2
- σc: uniaxial compressive strength
- σ11: negative of the x-coordinate of one vertex of cone (triangle in 2D) and negative of the y-coordinate of adjacend vertex in principal stress space

Note that the yield curve does not depend on the floe field. The `MohrsCone` object then has a `vertices` field that holds the coordiantes of the cone in principle stress space.


#### Simplification Settings
`SimplificationSettings` changes Subzero's floe simplification. The availible settings are:
- smooth_vertices_on, which turns on and off the process of smoothing floe vertices to decrease the total number of vertices
- max_vertices, the total number of verticies a floe can have before smoothing
- Δt_smooth, which sets the number of timesteps between floe smoothing
- tol, which is the tolerance in Douglas–Peucker's polygon simplification algorithm in meters.

Here is an example of creating your own simplification settings, using the default values:
```julia
simp_settings = SimplificationSettings(
    FT,
    min_floe_area = 1e6,
    smooth_vertices_on = true,
    max_vertices = 30,
    Δt_smooth = 20,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Ridging and Rafting Settings
`RidgeRaftSettings` changes Subzero's ridging and rafting process. The availible settings are:
- ridge_raft_on, which turns on and off the process of smoothing floe vertices to decrease the total number of vertices
- Δt  sets the number of timesteps between floe smoothing
- ridge_probability is the likelyhood (between 0-1) that two floes that meet the criterion actually ridge
- raft_probability is the likelyhood (between 0-1) that two floes that meet the criterion actually raft
- min_overlap_frac is the minimum overlap area fraction between a floe and another floe/domain for that floe to ridge or raft
- min_ridge_height is the minimum floe height to ridge with a floe/domain
- max_floe_ridge_height is the maximum floe height to ridge with another floe
- max_domain_ridge_height is the maximum floe height to ridge with a domain element
- max_floe_raft_height is the maximum floe height to raft with another floe
- max_domain_raft_height is the maximum floe height to raft with a domain element
- domain_gain_probability is the probalility that a floe that rafts with a domain element keeps all of its mass (0) or if that mass is removed and lost to the domain element (1).

Here is an example of creating your own simplification settings, using the default values:

```julia
ridgeraft_settings = RidgeRaftSettings(
  ridge_raft_on = false,
  Δt = 0,
  ridge_probability = 0.95,
  raft_probability = 0.95,
  min_overlap_frac = 0.01,
  min_ridge_height = 0.2,
  max_floe_ridge_height = 5.0,
  max_domain_ridge_height = 1.25,
  max_floe_raft_height = 0.25,
  max_domain_raft_height = 0.25,
  domain_gain_probability = 1.0,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Welding Settings
`WeldSettings` changes Subzero's welding process. The availible settings are:
- weld_on is a boolean flag for if welding should be turned on in the simulation
- Δts is a a list of multiples of timesteps during which welding code will run, welding will be run at multiples of all elements, each with domain split into corresponding Nx and Ny values
- Nxs is a list of number of x-directional bins to split the domain into at corresponding timesteps
- Nys is a list of number of x-directional bins to split the domain into at corresponding timesteps
- min_weld_area is the minimum area a weld can create for two floes to weld
- max_weld_area is the maximum area a weld can create for two floes to weld
- welding_coeff is a non-dimensional parameter, multiplied by ratio of overlap between two floes to original floe area to determin probability that a floe will merge. The larger this is, the more likely floes are to weld. Probability with 5% overlap is `welding_coeff * (0.05) > rand()`

Here is an example of creating your own welding settings, using the default values:

```julia
weld_settings = WeldSettings(
    weld_on = false,
    Δts = Vector{Int}(),
    Nxs = Vector{Int}(),
    Nys = Vector{Int}(),
    min_weld_area = 1e6,
    max_weld_area = 2e9,
    welding_coeff = 150,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

### Timesteps
You have the ability to set the simulation's timestep in seconds using `∆t` and set the total number of timsteps the simulation will run for, `n∆t`. The default is `∆t = 10` seconds and `n∆t = 7500` timesteps. 

### Output Writers
You can add four types of output writers, and as many of each type as you would like. The four types are as follows: `InitialStateOutputWriter`, `CheckpointOutputWriter`, `FloeOutputWriter`, and `GridOutputWriter`. When any of these objects are created, the file that they will write to is also created automatically. A brief desctiption of each is below:
#### InitialStateOutputWriter
The initial state output writer allows you to save the initial state of your simulation so that you can re-load it later. It saves the simulation object that you are currently creating to a JLD2 file. An `InitialStateOutputWriter` has two fields:
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `InitialStateOutputWriter` directly as a struct or with the following function call:
```julia
init_writer = InitialStateOutputWriter(
    dir = ".",
    filename = "initial_state.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
Note that these are the default values and you only need to pass in arguments that you wish to change from these values. Also note that all inputs are keyword arguments.

#### CheckpointOutputWriter
The checkpoint output writer allows you to save the state of the floes, ocean, and atmosphere at a specified part of the simulation into a JLD2 file. This will give you the ability to easily restart your simulation from the last-saved checkpoint if it were to fail for any reason, or if you wanted to continue a previous run from its stopping point. A `CheckpointOutputWriter` has three fields:
- `Δtout`, which specifies the umber of timesteps between checkpoint outputs starting from the first timestep
- `filename`, including the path to the file, to save the file to
-  `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `CheckpointOutputWriter` directly as a struct or with the following function call:
```julia
checkpointer = CheckpointOutputWriter(
    Δtout,
    dir = ".",
    filename = "checkpoint.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
Note that other than `∆tout`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### FloeOutputWriter
The floe output writer allows you to save floe values into a JLD2 file. This will give you the ability to easily analyze floe fields across timesteps. A `FloeOutputWriter` has four field:
- `Δtout`, which specifies the umber of timesteps between floe outputs starting from the first timestep
- `outputs`, which specifies which floe fields should be included
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `FloeOutputWriter` directly as a struct or with the following function call:
```julia
floewriter = FloeOutputWriter(
    Δtout;
    outputs = collect(fieldnames(Floe)),
    dir = ".",
    filename = "floes.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
The `outputs` field takes in a list of symbols corresponding to floe fields. For example, if you want the floe output writer to output the floes centroid and coordinates then `outputs = [:centroid, :coords]`. If you want all floe fields then you can simply omit the outputs field all together and all floe fields will be output. Note that other than `∆tout`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

Note that if you have Periodic calls, and thus ghost floes in your simulation, these will also be saved by the `FloeOutputWriter`. If you want to exclude these floes from your analysis or when otherwise using the `FloeOutputWriter` output, you can do so by only including floes with a `ghost_id = 0`.

#### GridOutputWriter
The grid output writer allows you to floe values averaged onto a course grid to a NetCDF file. This will give you the ability to easily analyze floe characteristics on a grid. A `GridOutputWriter` has eight field:
- `outputs`, which specifies which floe fields should be included
- `Δtout`, which specifies the umber of timesteps between floe outputs starting from the first timestep
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.
- `xg`, the grid lines in the x-direction of the grid that you would like the calculations done over (doesn't have to be the same as simulation grid)
- `yg`, the grid lines in the y-direction of the grid that you would like the calculations done over (doesn't have to be the same as simulation grid)
- `data`, three-dimensional array that holds grid averaged data prior to writing to file
- `average`, boolean that specifies if the gridded data should be averaged over each timestep between writing to file, or if it should just be calculated prior to writing for that singular timestep (NOT IMPLEMENTED YET).

You can create an `GridOutputWriter` directly as a struct or with the following function call:
```julia
gridwriter = GridOutputWriter(
    FT,
    Δtout,
    grid,
    dims;
    outputs = collect(get_known_grid_outputs()),
    dir = ".",
    filename = "gridded_data.nc",
    overwrite = false,
    average = false,
)
```
The `outputs` field takes in a list of symbols. To see all possible outputs, call the `get_known_grid_outputs()` function. For example, if you want the grid output writer to output the floe masses and areas averaged on the grid then `outputs = [:mass_grid, :area_grid]`. If you want all possible fields then you can simply admit the outputs field altogether and all grid fields will be output. The `grid` field is the simulation grid, and then `dims` field specifies the dimensions of the grid you would like the output calculated on.

Note that other than `∆tout`, `grid`, and `dims`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### OutputWriters
Once you have created all of the types of output writers you need, you must combine them into one `OutputWriters` object that will be a simulation field. 
The `OutputWriters` struct has four fields: `initialwriters`, `floewriters`, `gridwriters`, and `checkpointwriters`. For each you can supply a StructArray of the specified type of writer. If you do not have any writers of a given type, don't provide any. Below I create two example `OutputWriters` objects, assuming that I have already created the following outputwriters: `initwriter1`, `checkpointer1`, `floewriter1`, `floewriter2`, and `gridwriter1`.

```julia
using StructArrays, Subzero

outputwriters1 = OutputWriters(initwriter1, floewriter1)

outputwriters2 = OutputWriters(
   initwriter1,
   checkpointer1,
   floewriter1,
   floewriter2,
   gridwriter1,
)
```
Here you can see that you can choose which values to supply and that you can supply more than one of each type if desired. You might want to do this is you want different outputs at different timeframes. 

### Reproducibility
The simulations are currently completly reproducible when run single-threaded. The simulation takes a `rng` argument, and if it is provided with a seeded random number generator, the two simulations with the same set of starting floes will produce the exact same results. Note that the same set of floes can be reproduced by providing a seeded random number generator to the floe creation functions. 

This has not yet been achieved for multi-threaded runs. 

### Creating the Simulation
Once you have created all of the above objects, you can combine them to create a `Simulation`. A simulation has quite a few fields, all of which are talked about above in more detail. Since there are so many fields, they are keyword defined, so you must provide a keyword when creating the struct. The only necessary argument is the `model` as everything else has a default value. However, a table of optional elements, and their default values, is as follows:

| Keyword           |    Default Value        | What is it?                                               |
| ----------------- | ------------------------| --------------------------------------------------------- |
| consts            | Constants()             | Physical parameters used in the simulation                |
| rng               | Xoshiro()               | Random number generator - can seed for reproducibility    |
| verbose           | false                   | Flag for printing timesteps and updates during simulation |
| name              | "sim"                   | Name of simulation                                        |
| ∆t                | 10                      | Length of simulation timestep in seconds                  |
| n∆t               | 7500                    | Number of timesteps run in simulation                     |
| floe_settings     | FloeSettings()          | Settings for making new floes    during the simulation    |
| coupling_settings | CouplingSettings()      | Settings for coupling during the simulation               |
| collision_settings| CollisionSettings()     | Settings for collisions during the simulation             |
| fracture_settings | FractureSettings()      | Settings for fractures during the simulation              |
| simp_settings     | SimplificationSettings()| Settings for floe simplification during the simulation    |
| ridgeraft_settings| RidgeRaftSettings()     | Settings for ridging and rafting during the simulation    |
| weld_settings     | WeldSettings()          | Settings for welding during the simulation                |
| writers           | OutputWriters()         | Lists of output writers to be written to during the run   |

Here is an example of how to create a simulation, providing some of these optional inputs given that we have already created the following objects: `my_model`, `my_consts`, `my_coupling_settings`, `my_simp_settings`, and `my_writers`.

```julia
my_simulation = Simulation(
    model = my_model,
    consts = my_consts,
    Δt = 5,  # Timestep of 5 seconds
    coupling_settings = my_coupling_settings,
    simp_settings = my_simp_settings,
    writers = my_writers,
)
```

### Running the simulation

You can now use the `run!` function to run the simulation:
```
run!(my_simulation)
```

If you wish to couple to Oceananigans, you will need to run each model timestep by timestep and pass the needed fields back and forth. You can run a single timestep of the simulation using the `timestep_sim!` function. This also needs the current timestep the simulation is on (`tstep`) as an argument.

```
timestep_sim!(
   my_simulation,
   tstep,
)
```

Note that we are working on a more elegant solution to coupling with Oceananigans and CliMA and this page will be updated once that is in place. 

If you run your simulation in multiple parts and need to re-start your simulation from files, the `restart!` function will be a good place to start. However, note that it is quite simple and users may need to write their own restart function if they want any complex behavior. 

The provided `restart!` function takes in the output file from both an `InitialStateOutputWriter` and a `CheckpointOutputWriter` to restart the simulation. In addition to providing these two files, the user must also provide the number of timesteps to run the next part of the simulation for (`new_nΔt`) and new output writers. The user also has an option to specify a non-zero starting timestep for the simulation using the keyword argument `start_tstep`.

```
restart!(
   initial_state_fn,
   checkpointer_fn,
   new_nΔt,
   new_output_writers;
   start_tstep = 0,
)
```

### Plotting

If your simulation has both a `FloeOutputWriter` and an `InitialStateOutputWriter`, you can use the built in plotting function to make an MP4 file with each frame as a timestep saved by the `FloeOutputWriter`. You do this as follows:

```
plot_sim(
    floe_output_writer_file_path,
    initial_state_output_writer_file_path,
    Δt,
    output_file_path,
)
```

where `floe_output_writer_file_path` is the .jl file saved by the `FloeOutputWriter`, `initial_state_output_writer_file_path` is the .jl file saved by the `InitialStateOutputWriter`, and `output_file_path` is the file path and name you want your .mp4 file saved as. `Δt` is the model timestep.
