"""Floes for comparing collision simulations to MATLAB - see matching floes in example_sims.m"""
FT = Float64
# (1) Two floes hitting straight on at close range - Gain ~ 1% energy --> matches sim (1) in matlab file
floe1_poly = LG.Polygon([[[1.25e4, 8e4], [1.25e4, 6e4], [1.75e4, 6e4], 
                    [1.75e4, 8e4], [1.25e4, 8e4]]])
floe1 = Floe(floe1_poly, hmean, Δh, u = 0.0)
floe2_poly = LG.Polygon([[[0.7e4, 8e4], [0.7e4, 6e4], [1.2e4, 6e4], 
                    [1.2e4, 8e4], [0.7e4, 8e4]]])
floe2 = Floe(floe2_poly, hmean, Δh, u = 1.0)
floe_arr = StructArray([floe1 floe2])

# (2) Two floes hitting - one really big and one small -  they do not sink into each other like the edges of the simulation
floe1_poly = LG.Polygon([[[1.25e4, 9e4], [1.25e4, -9e4], [3.75e4, -9e4], 
                    [3.75e4, 9e4], [1.25e4, 9e4]]])
floe1 = Floe(floe1_poly, hmean, Δh, u = 0.0)
floe2_poly = LG.Polygon([[[0.7e4, 8e4], [0.7e4, 6e4], [1.2e4, 6e4], 
                    [1.2e4, 8e4], [0.7e4, 8e4]]])
floe2 = Floe(floe2_poly, hmean, Δh, u = 1.0)
floe_arr = StructArray([floe1 floe2])

# Floe hitting the boundary - zero ocean, zero wind, needs friction
floe1_poly = LG.Polygon([[[9.25e4, 7e4], [9.25e4, 5e4], [9.75e4, 5e4], 
                    [9.75e4, 7e4], [9.25e4, 7e4]]])
floe1 = Floe(floe1_poly, hmean, Δh, u = 1.0)
floe_arr = StructArray([floe1])

# Floe hitting topography
topo = TopographyElement(; coords = [[[0.5e4, 5e4], [0.5e4, 7e4], [1e4, 7e4], [1e4, 5e4], [0.5e4, 5e4]]])
topo_arr = StructVector([topo for i in 1:1])

domain = Subzero.Domain(; north = nboundary, south = sboundary, east = eboundary, west = wboundary, topography = topo_arr)

floe1_poly = LG.Polygon([[[-1.75e4, 5e4], [-1.75e4, 7e4], [-1.25e4, 7e4], 
                    [-1.25e4, 5e4], [-1.75e4, 5e4]]])
floe1 = Floe(floe1_poly, hmean, Δh, u = 1.0)
floe_arr = StructArray([floe1])