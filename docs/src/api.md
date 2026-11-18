```@meta
CurrentModule = Subzero
```

# Full Subzero API documentation

## User Methods

### Grids

```@docs
AbstractRectilinearGrid
RegRectilinearGrid
```

### Directions
```@docs
AbstractDirection
North
South
East
West
```
### Boundaries
```@docs
AbstractBoundary
OpenBoundary
PeriodicBoundary
CollisionBoundary
MovingBoundary
```
### Topography
```@docs
TopographyElement
initialize_topography_field
```
### Domain
```@docs
Domain
```

### Ocean
```@docs
Ocean
```

### Atmosphere
```@docs
Atmos
```

### Floes
```@docs
Floe
AbstractFloeFieldGenerator
CoordinateListFieldGenerator
VoronoiTesselationFieldGenerator
initialize_floe_field
```

### Model
```@docs
Model
```

### Constants
```@docs
Constants
```

### Physical Process Settings
```@docs
FloeSettings
CouplingSettings
CollisionSettings
FractureSettings
AbstractFractureCriteria
NoFracture
HiblerYieldCurve
MohrsCone
SimplificationSettings
RidgeRaftSettings
WeldSettings
```

### Output Writers
```@docs
AbstractOutputWriter
InitialStateOutputWriter
CheckpointOutputWriter
GridOutputWriter
FloeOutputWriter
OutputWriters
SubzeroLogger
```

### Simulations
```@docs
Simulation
run!
restart!
```

## Developer Methods

### Simulation Methods
```@docs
timestep_sim!
```
### Output Writer Methods
```@docs
write_data!
write_init_state_data!
write_checkpoint_data!
write_floe_data!
write_grid_data!
calc_eulerian_data!
```
### Developer-Used Types
```@docs
CellFloes
CellStresses
TopographyField
```

### Collision Methods
```@docs
timestep_collisions!
calc_elastic_forces
calc_normal_force
calc_friction_forces
add_interactions!
floe_floe_interaction!
floe_domain_interaction!
floe_domain_element_interaction!
update_boundaries!
add_ghosts!
```

### Coupling Methods
```@docs
timestep_coupling!
calc_one_way_coupling!
calc_two_way_coupling!
calc_atmosphere_forcing
calc_ocean_forcing!
floe_to_grid_info!
```

### Fracture Methods
```@docs
fracture_floes!
split_floe
deform_floe!
determine_fractures
```

### Ridging and Rafting Methods
```@docs
timestep_ridging_rafting!
floe_floe_raft!
floe_domain_raft!
floe_floe_ridge!
floe_domain_ridge!
```

### Simplify Methods
```@docs
simplify_floes!
smooth_floes!
fuse_floes!
remove_floes!
```

### Welding Methods
```@docs
timestep_welding!
```

### Floe Update Methods
```@docs
timestep_floe_properties!
calc_stress!
calc_strain!
```
