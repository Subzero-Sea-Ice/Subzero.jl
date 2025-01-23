```@meta
CurrentModule = Subzero
```

# Full Subzero API documentation

!!! warning
    This page is still very much WIP!

## Grids

```@docs
AbstractRectilinearGrid
RegRectilinearGrid
```

## Directions
```@docs
AbstractDirection
North
South
East
West
```
## Boundaries
```@docs
AbstractBoundary
OpenBoundary
PeriodicBoundary
CollisionBoundary
MovingBoundary
```
## Topography
```@docs
TopographyElement
initialize_topography_field
```
## Domain
```@docs
Domain
```

## Ocean
```@docs
Ocean
```

## Atmosphere
```@docs
Atmos
```

### Model
```@docs
Model
```

## Developer-Used Types
```@docs
CellFloes
CellStresses
TopographyField
```
