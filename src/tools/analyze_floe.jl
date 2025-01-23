# Needs to stop using Plots if used again

f = jldopen("output/voronoi/fresh_start_2246.jld2", "r")

u_arr = Vector{Float64}()
for k in range(0, 1470, step = 10)
    push!(u_arr, f[string("u/", k)][1])
end

collision_arr = Vector{Float64}()
for k in range(0, 1470, step = 10)
    push!(collision_arr, f[string("collision_force/", k)][1])
end

mass_arr = Vector{Float64}()
for k in range(2246, 2360)
    push!(mass_arr, f[string("mass/", k)][1])
end

p_dudt_arr = Vector{Float64}()
for k in range(2246, 2360)
    push!(p_dudt_arr, f[string("p_dudt/", k)][1])
end

x_cforce = first.(first.(collision_arr))
fxOA_arr = first.(fxOA_arr)
total_force = fxOA_arr .+ x_cforce
total_force_per_mass = total_force ./ mass_arr
frac = p_dudt_arr ./ total_force_per_mass

plot(range(2246, 2360), u_arr)
savefig("output/voronoi/u.png")

plot(range(2246, 2360), x_cforce)
savefig("output/voronoi/collision_x.png")

plot(range(2246, 2360), first.(first.(centroid_arr)))
savefig("output/voronoi/centroid_x.png")

plot(range(2246, 2360), fxOA_arr)
savefig("output/voronoi/fxOA.png")

plot(range(2246, 2360), p_dudt_arr)
savefig("output/voronoi/p_dudt.png")

plot(range(2246, 2360), fxOA_arr .+ x_cforce)
savefig("output/voronoi/total_force.png")

plot(range(2246, 2360), frac)
savefig("output/voronoi/frac.png")

close(f)

#  plotting voronoi w/ labels
scatter(xp, yp, markersize = 3, label = "centroids")
annotate!([(xp[n] + 0.2, yp[n] + 0.03, Plots.text(id[n], 4)) for n in 1:159])
xs, ys = first.(c[1]), last.(c[1])
plot!(xs, ys, seriestype = [:shape])
savefig("output/sim/labeled_voronoi.png")