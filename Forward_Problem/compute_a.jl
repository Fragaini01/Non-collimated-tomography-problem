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
using PyCall
using CSV, DataFrames
using NearestNeighbors
using JLD2


interpolate = pyimport("scipy.interpolate")

include("Functions.jl")

# Variables
verbose = false
plot = false
test = true

filepath_equil = "equilibrium/ITER/test/80MW_equilibrium.geqdsk"

if verbose
    println("Reading equilibrium file: ", filepath_equil)
end

M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction

if verbose
    println("Equilibrium file read.")
end

R_hfs = minimum(wall.r) # R-coord of high-field side wall
R_lfs = maximum(wall.r) # R-coord of low-field side wall
println("R_lfs: ", R_lfs)
wall_phi = collect(0:1:359).*(2*pi/180.0) # Toroidal angle
topview_R_hfs_x = (R_hfs).*cos.(wall_phi)
topview_R_hfs_y = (R_hfs).*sin.(wall_phi)
topview_R_lfs_x = (R_lfs).*cos.(wall_phi)
topview_R_lfs_y = (R_lfs).*sin.(wall_phi)


first_wall_size = 0.15
second_wall_size = 0.40
R = 6.5
cross_section_1 = 11.1
cross_section_2 = 19.4
radii, zs = regrid_wall(wall , .001)
regrgridded_wall = Equilibrium.Boundary(radii, zs)

plasma_grid = plasma_gridder(.1, 1, .1, radii, zs)
grid_henrik = plasma_gridder_henrik(100, 200, wall)

detectors_file = "../my_code/detectors.csv"
detectors_df = CSV.read(detectors_file, DataFrame; header=true)

detectors = [collect(row) for row in eachrow(detectors_df)]

println("Number of detectors: ", length(detectors))

first_slice = plasma_grid[1:Int(length(plasma_grid) / 359)]
slice_size = Int(length(first_slice))


# A calculator 
a = zeros(Float64, length(first_slice), length(detectors))
time_a = 0.0
total_start = time()
points_in_the_wall = []


for j in 1:length(detectors)
    println("Detector $j: $(round(detectors[j][1], digits=2)) , $(round(detectors[j][3], digits=2))")
    # Leggi i voxel visti dal detector j dal file output/detector_j
    filename = "../my_code/output/detector_$(j).vc"
    global voxels_seen = []
    if isfile(filename)
        df = CSV.read(filename, DataFrame; header=false,)
        global voxels_seen = [[df[i,12], df[i,13], df[i,6]] for i in 1:nrow(df)]
        #global voxels_solid = [df[i , end] for i in 1:nrow(df)]
    else
        println("File not found: ", filename)
        abort()
    end
    if verbose
        println("Number of voxels seen by detector $(j): ", length(voxels_seen))
    end
    t_start = time()
    global mask = mask_grid_interpolate(voxels_seen, grid_henrik, plasma_grid , regrgridded_wall)
    t_end = time()
    global time_a += (t_end - t_start)
    if verbose
        println("Mask grid interpolation time: ", t_end - t_start)
    end

    first_guess = [rand(1:length(radii)), rand(1:length(radii))]
    for point in 1:length(first_slice)
        total_attenuation = 0 
        for angle in 0:358
            i = point + slice_size*angle
            if mask[i] == 1
                appo = a_ij(plasma_grid[i], detectors[j], radii, zs, first_wall_size, second_wall_size, cross_section_1, cross_section_2 , first_guess)
                total_attenuation += appo[1]
                first_guess = [appo[2], appo[3]] 
            end
        end
        a[point, j] = total_attenuation 
    end
    if verbose
        println("Total time for detector $(j): ", time() - t_start)
    end
end

# Varius tests to ensure the matrix a is correct
if test 
    println("Shape of a: ", size(a))
    println("Number of detectors: ", size(a, 2))
    println("Number of points in first slice: ", size(a, 1))
    println("Columns with all zeros: ", findall(j -> all(==(0.0), a[:,j]), 1:size(a,2)))

    println("min(a) = ", minimum(vec(a)), ", max(a) = ", maximum(vec(a)))
    println("Number of zeros in a: ", count(==(0.0), a))
    println("Contains NaN? ", any(isnan, a), " | Contains Inf? ", any(isinf, a))
    tol = 1e-10
    rows_almost_zero = findall(i -> all(x -> abs(x) < tol, a[i,:]), 1:size(a,1))
    println("symmetric and positive $(issymmetric(a'a) && all(eigvals(a'a) .> 0))")
end

total_end = time()

println("Time to compute a: ", total_end - total_start)
println("mask_grid_interpolate execution time: ", time_a)

a_to_W(a , plasma_grid, detectors , "../my_code/output/a")

### Plotting
if plot
    psi_mag, psi_bdry = psi_limits(M)
    wall_dR = maximum(wall.r)-minimum(wall.r)
    plot_font = "Computer Modern"
    Plots.default(fontfamily=plot_font)

    # Loop over all detectors and create/save a plot for each
    nR = length(unique([p[1] for p in first_slice]))
    nZ = length(unique([p[3] for p in first_slice]))
    clims_common = (0.0, maximum(a))
    color = cgrad([:white, :blue, :red], [0.0, 0.5, 1.0])
    colorbar = true
    for det_idx in 1:length(detectors)
        a_to_plot = reshape(a[:,det_idx], nR, nZ)
        plt_crs = Plots.plot(aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]", xlims=(3, 11))
        plt_crs = Plots.scatter!(
            plt_crs,
            [p[1] for p in first_slice],         # R coordinate
            [p[3] for p in first_slice],         # z coordinate
            marker_z = vec(a_to_plot),           # color by a_to_plot values
            color = color,                       # diverging colormap
            clims = (0.0 , maximum(a_to_plot)),                # fixed color limits for all detectors
            markerstrokecolor = :transparent, 
            markerstrokewidth = 0,               # remove marker border
            #seriesalpha = abs.(vec(a_to_plot)) ./ maximum(abs, a_to_plot), # transparency proportional to |a_to_plot|            
            colorbar = colorbar,
            label = "Grid points",
            ms = 5, 
        )
        plt_crs = Plots.plot!(plt_crs, radii, zs, label="Wall spline", linecolor=:black, linewidth=1)
        plt_crs = Plots.scatter!(plt_crs, wall.r, wall.z, label="Tokamak first wall", linewidth=2.5, color=:black)
        # Plot only the current detector position on the cross-sectional view
        det_R = detectors[det_idx][1]
        det_z = detectors[det_idx][3]
        plt_crs = Plots.scatter!(
            plt_crs,
            [det_R],
            [det_z],
            marker=:star5,
            color=:red,
            ms=8,
            label="Detector $det_idx"
        )
        plt_crs = Plots.plot!(plt_crs, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16)
        plt_crs = Plots.plot!(plt_crs, legend=:bottomright, legendfontsize=13)
        plt_crs = Plots.plot!(plt_crs, title="Detector $det_idx", titlefontsize=14)
        plt_crs = Plots.plot!(plt_crs, dpi=200)
        plt_crs = Plots.plot!(plt_crs, xlabel="R [m]", ylabel="z [m]")
        plt_crs = Plots.plot!(plt_crs, legend=false)
        println("Saving plot for detector $det_idx...")
        savefig(plt_crs, "../my_code/output/plots/a_detector_$(det_idx).pdf")
    end

end