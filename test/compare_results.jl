using Makie, CairoMakie, NCDatasets, MAT
mat_path = "src/SubZero_MATLAB/Floes/"
mat_files = readdir(mat_path)
u_mat = Float64[]
v_mat = Float64[]
for file in mat_files
    mat_data = matread(joinpath(mat_path, file))
    push!(u_mat, mat_data["Floe"]["Ui"])
    push!(v_mat, mat_data["Floe"]["Vi"])
end

jul_path = "output/sim/f.nc"
jul_data = NCDataset(jul_path)
u_jul = jul_data["u"][:, 1]
v_jul= jul_data["v"][:, 1]
NCDatasets.close(jul_data)

u0 = 1.0
v0 = 0.0
ρi = 920
ρo = 1027
Cd = 3e-3
h = 0.25
τu = (ρi*h)/(ρo*Cd*abs(u0))
τv = (ρi*h)/(ρo*Cd*abs(v0))
Δt = 10  # seconds
t = collect(0:50:3000) .* Δt
theoretical_u = u0 * (1 .- (1 ./(t./τu.+1)))
theoretical_v = v0 * (1 .- (1 ./(t./τv.+1)))
t_min = t./60


scatter(t_min, theoretical_u, title = "1m/s Zonal and 0.0m/s Meridional Flow:\nU Velocity Theoretical vs MATLAB vs Julia", label = "Theoretical", legend=:bottomright, xlabel = "Time [min]", ylabel = "Velocity [m/s]")
scatter!(t_min, u_jul, label = "Julia")
scatter!(t_min[2:end], u_mat, label = "Matlab")
savefig("src/jmt_ucompare_north_west")

scatter(t_min, theoretical_v, title = "-0.5m/s Zonal and 0.5m/s Meridional Flow:\nV Velocity Theoretical vs MATLAB vs Julia", label = "Theoretical", legend=:bottomright, xlabel = "Time [min]", ylabel = "Velocity [m/s]")
scatter!(t_min, v_jul, label = "Julia")
scatter!(t_min[2:end], v_mat, label = "Matlab")
savefig("src/jmt_vcompare_north_west")