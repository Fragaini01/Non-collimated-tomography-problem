using Distributed 
folderpath_OWCF = "/home/francesco/Documenti/Master/Blok4/PUK/OWCF/"

cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using Dates
using FileIO
using LinearAlgebra
using Equilibrium
using Plots #; gr()
using Random
using Dierckx
using Dates
using CSV, DataFrames
using NearestNeighbors
using JLD2

include(folderpath_OWCF*"misc/temp_n_dens.jl")

include("Functions.jl")

# Variables
verbose = false
plot = false

filepath_equil = "equilibrium/ITER/test/80MW_equilibrium.geqdsk"

if verbose
    println("Reading equilibrium file: ", filepath_equil)
end

M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction

if verbose
    println("Equilibrium file read.")
end

dR = 0.1
dz = 0.1
radii, zs = regrid_wall(wall , .01)
plasma_grid = plasma_gridder(dR, 1, dz, radii, zs)

center = [magnetic_axis(M)[1], magnetic_axis(M)[2]]
slice_size = Int(length(plasma_grid) / 359)
println("Slice size: ", slice_size)

first_slice = plasma_grid[1:slice_size]
# In Julia, the equivalent of assert is the @assert macro, which is built-in and does not require any import.
# Replace 'assert' with '@assert' as follows:

@assert length(first_slice) == slice_size "First slice size does not match expected size."

#Assuming all the plasma's slices are the same, and the only distance is on R-Z plane
neutros = zeros(Float64, length(first_slice))
println("Number of neutros: ", length(neutros))

idx = 1
for voxel in first_slice
    if in_boundary(wall , voxel[1], voxel[3])
        psi_rz = M(voxel[1], voxel[3])
        psi_on_axis, psi_at_bdry = psi_limits(M)
        rho_pol_rz = sqrt(max(0.0,(psi_rz-psi_on_axis)/(psi_at_bdry-psi_on_axis)))

        T_i = getAnalyticalTemp(10.0, rho_pol_rz)

        n_i = getAnalyticalDens(1.0e20, rho_pol_rz)
        C_1 = 1.17*10^(-9);C_2 = 1.51*10^(-2);C_3 = 7.52*10^(-2);C_4 = 4.61*10^(-3);C_5 = 1.35*10^(-2);C_6 = -1.07*10^(-4);C_7 = 1.37*10^(-5);
        B_G = 34.4

        theta = T_i/ (1 .- T_i* (C_2 .+ T_i* (C_4 .+ T_i* C_6)) ./ (1 .+ T_i* (C_3 .+ T_i* (C_5 .+ T_i* C_7))));

        xi = (B_G^2 ./ (4 .* theta)).^(1/3)

        reac = C_1*theta*(xi/(1.125*10^6*T_i^3))^(1/2)*exp(-3*xi);
        
        value = n_i * T_i
        if value >  (3/3.7)*10^21
            neutros[idx] = n_i * reac
        end
    end
    global idx += 1
end

# Prepara il plot di sfondo (parete)
psi_mag, psi_bdry = psi_limits(M)
wall_dR = maximum(wall.r)-minimum(wall.r)
plot_font = "Computer Modern"
Plots.default(fontfamily=plot_font)
vmax = maximum(vec(neutros))

radii = [p[1] for p in first_slice]
zs = [p[3] for p in first_slice]

# Prima scatter: true values
scatter(
    vec(radii), 
    vec(zs),
    marker_z = vec(neutros),
    color = :balance,
    clims = (-vmax, vmax),
    ms = 6,
    markerstrokecolor = :transparent,
    title = "True Intensities",
    xlabel = "R",
    ylabel = "Z",
    colorbar = true,
    aspect_ratio = :equal,
    label = "", 
    xlims = (3, 9)
)


println("Number of neutros: ", length(neutros))

for n in [20,50,99]
    compute_measurements(neutros, "../my_code/output/a" , "../my_code/output/measurements" , n)
end