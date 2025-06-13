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
using Measures
plot_seen = false
plot_a = false
plot_all_detector = false
plot_grid = false
plot_suma = true

verbose = true

interpolate = pyimport("scipy.interpolate")

include("Functions.jl")


filepath_equil = "equilibrium/ITER/test/80MW_equilibrium.geqdsk"

if verbose
    println("Reading equilibrium file: ", filepath_equil)
end

M, wall = read_geqdsk(filepath_equil,clockwise_phi=false) # Assume counter-clockwise phi-direction

if verbose
    println("Equilibrium file read.")
end
# Variables
first_wall_size = 0.15
second_wall_size = 0.40
R = 6.5
cross_section_1 = 11.1
cross_section_2 = 19.4
radii, zs = regrid_wall(wall , .001)
regrgridded_wall = Equilibrium.Boundary(radii, zs)

plasma_grid = plasma_gridder(.1, 1, .1, radii, zs)
grid_henrik = plasma_gridder_henrik(100, 200, wall)

detectors_file = "../my_code/detectors(1).csv"
detectors_df = CSV.read(detectors_file, DataFrame; header=true)

detectors = [collect(row) for row in eachrow(detectors_df)]

println("Number of detectors: ", length(detectors))

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

first_slice = plasma_grid[1:Int(length(plasma_grid) / 359)]
slice_size = Int(length(first_slice))

# read W from file
in_file = "/home/francesco/Documenti/Master/Blok4/PUK/my_code/output/a.jld2"

if isfile(in_file )
    @load (in_file ) W Ed_array Ed_array_units D1_array D1_array_units D2_array D2_array_units
else
    println("File not found: ", in_file )
    return
end
if length(W) == 0
    println("W is empty, cannot compute measurements.")
    return
end
W = reshape(W, length(Ed_array), length(D1_array) * length(D2_array))

a = W'


if plot_seen
    for j in 1:length(detectors)
        println("Detector $j seen points")
        detector_location = detectors[j]

        filename = "../my_code/output/detector_$(j).vc"

        if isfile(filename)
            df = CSV.read(filename, DataFrame; header=false,)

            voxels_seen = [[df[i,12], df[i,13]] for i in 1:nrow(df)]
            voxels_seen = unique(voxels_seen)
        else
            println("File not found: ", filename)
            abort()
        end


        plt_top = Plots.plot(topview_R_hfs_x,topview_R_hfs_y,linewidth=2.5,color=:black,label="")
        plt_top = Plots.plot!(topview_R_lfs_x,topview_R_lfs_y,linewidth=2.5,color=:black,label="")
        plt_top = Plots.plot!(aspect_ratio=:equal,xlabel="x [m]",ylabel="y [m]")
        plt_top = Plots.plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)
        plt_top = scatter!(
            [p[1] * cos(p[2]) for p in voxels_seen],
            [p[1] * sin(p[2]) for p in voxels_seen],
            markersize=2,
            label="",
            color=:red,
            xlabel="x [m]",
            ylabel="y [m]",
            legend=:bottomright,
            aspect_ratio=:equal,
            title="Voxels seen by Detector $j",
            markerstrokecolor = :transparent, 
            markerstrokewidth = 0,   
        )
        plt_top = Plots.scatter!([detector_location[1]],[detector_location[2]], label="Detector location", markershape=:star, markercolor=:purple, markerstrokewidth=2)

        savefig(plt_top, "../my_code/output/plots/detector_$(j).pdf")
    end
end




if plot_a
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
        plt_crs = Plots.plot(
            aspect_ratio = :equal,
            xlabel = "R [m]",
            ylabel = "z [m]",
            xlims = (3, 12),
            size = (900, 700),           # aumenta la dimensione
            left_margin = 5mm,           # margini piccoli
            right_margin = 15mm,         # più spazio a destra per la colorbar
            top_margin = 5mm,
            bottom_margin = 5mm,
        )

        plt_crs = Plots.scatter!(
            plt_crs,
            [p[1] for p in first_slice],         # R coordinate
            [p[3] for p in first_slice],         # z coordinate
            marker_z = vec(a_to_plot),           # color by a_to_plot values
            color = color,                       # diverging colormap
            clims = (0.0 , maximum(a_to_plot)),                # fixed color limits for all detectors
            markerstrokecolor = :transparent, 
            markerstrokewidth = 0,               # remove marker border
            colorbar = colorbar,
            label = "Grid points",
            colorbar_padding = 40,
            ms = 5, 
            colorbar_title = "a",
            colorbar_titlefont = 16,
            colorbar_tickfont = 16,
            colorbar_orientation = :vertical,
            colorbar_position = :right, 
        )

        plt_crs = Plots.plot!(plt_crs, radii, zs, label="Wall spline", linecolor=:black, linewidth=2)
        plt_crs = Plots.scatter!(plt_crs, wall.r, wall.z, label="Tokamak first wall", linewidth=4, color=:black)
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
        plt_crs = Plots.plot!(plt_crs, grid = true)
        plt_crs = Plots.plot!(plt_crs, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16)
        plt_crs = Plots.plot!(plt_crs, legend=:bottomright, legendfontsize=13)
        plt_crs = Plots.plot!(plt_crs, title="Detector $det_idx", titlefontsize=18)
        plt_crs = Plots.plot!(plt_crs, dpi=200)
        plt_crs = Plots.plot!(plt_crs, xlabel="R [m]", ylabel="z [m]")
        plt_crs = Plots.plot!(plt_crs, legend=false)
        println("Saving plot for detector $det_idx...")
        savefig(plt_crs, "../my_code/output/plots/a_detector_$(det_idx).pdf")
    end
end


if plot_all_detector
    # Plot the reactor wall and all detectors in the R-z plane (cross-sectional view)
    plt_detectors_rz = Plots.plot(
        radii, zs,
        linewidth=2.5, color=:black, label="",
        aspect_ratio=:equal, xlabel="R [m]", ylabel="z [m]",
        title="Reactor wall with detectors (R-z plane)",
        xlims = (3, 12),

    )
    Plots.scatter!(plt_detectors_rz, wall.r, wall.z, label="", color=:black)

    # Add all detectors as stars in R-z
    Plots.scatter!(
        plt_detectors_rz,
        [d[1] for d in detectors],
        [d[3] for d in detectors],
        markershape=:star5, markercolor=:purple, ms=8, label="Detectors"
    )

    Plots.plot!(plt_detectors_rz, legend=:bottomright, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16)
    savefig(plt_detectors_rz, "../my_code/output/plots/all_detectors_rz(1).pdf")
    println("Saved plot of all detectors in R-z plane.")
end


if plot_suma
    # Plot della somma di tutte le righe di a
    sum_a = sum(a, dims=2)
    nR = length(unique([p[1] for p in first_slice]))
    nZ = length(unique([p[3] for p in first_slice]))
    sum_a_reshaped = reshape(sum_a, nR, nZ)

    plt_sum = Plots.plot(
        aspect_ratio = :equal,
        xlabel = "R [m]",
        ylabel = "z [m]",
        xlims = (3, 12),
        size = (900, 700),
        left_margin = 5mm,
        right_margin = 15mm,
        top_margin = 5mm,
        bottom_margin = 5mm,
        title = "Sum over all rows of a"
    )

    plt_sum = Plots.scatter!(
        plt_sum,
        [p[1] for p in first_slice],
        [p[3] for p in first_slice],
        marker_z = vec(sum_a_reshaped),
        color = cgrad([:white, :blue, :red], [0.0, 0.5, 1.0]),
        clims = (0.0, maximum(sum_a_reshaped)),
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        colorbar = true,
        ms = 5,
        colorbar_title = "Σ a",
        colorbar_titlefont = 16,
        colorbar_tickfont = 16,
        colorbar_orientation = :vertical,
        colorbar_position = :right,
    )

    plt_sum = Plots.plot!(plt_sum, radii, zs, label="Wall spline", linecolor=:black, linewidth=2)
    plt_sum = Plots.scatter!(plt_sum, wall.r, wall.z, label="Tokamak first wall", linewidth=4, color=:black)
    plt_sum = Plots.plot!(plt_sum, grid = true)
    plt_sum = Plots.plot!(plt_sum, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16)
    plt_sum = Plots.plot!(plt_sum, legend=false)
    plt_sum = Plots.plot!(plt_sum, dpi=200)
    plt_sum = Plots.plot!(plt_sum, xlabel="R [m]", ylabel="z [m]")

    savefig(plt_sum, "../my_code/output/plots/sum_a(1).pdf")
    println("Saved plot of sum of all rows of a in R-z plane.")
end


if plot_grid
    plasma_grid = plasma_gridder(.1, 1, .1, radii, zs)
    grid_henrik = plasma_gridder_henrik(100, 200, wall)

    radii_k2, zs_k2 = regrid_wall(wall, .001 , 2)
    radii_k3, zs_k3 = regrid_wall(wall, .001 , 3)
    radii_k4, zs_k4 = regrid_wall(wall, .001 , 4)

    ## Plotting
    plot_font = "Computer Modern"
    Plots.default(fontfamily=plot_font)
    plt = plot(layout = (2,2), size=(900,900), xlabel="R [m]", ylabel="z [m]", legend=false,     aspect_ratio = :equal,
    )

    plot!(plt[1], radii, zs, label="order=1", title="Spline order 1", color=:blue)
    scatter!(plt[1], wall.r, wall.z, label="Wall points", color=:black)

    plot!(plt[2], radii_k2, zs_k2, label="order=2", title="Spline order 2", color=:red)
    scatter!(plt[2], wall.r, wall.z, label="Wall points", color=:black)

    plot!(plt[3], radii_k3, zs_k3, label="order=3", title="Spline order 3", color=:green)
    scatter!(plt[3], wall.r, wall.z, label="Wall points", color=:black)

    plot!(plt[4], radii_k4, zs_k4, label="order=4", title="Spline order 4", color=:purple)
    scatter!(plt[4], wall.r, wall.z, label="Wall points", color=:black)

    savefig(plt, "../my_code/output/plots/wall_comparison.pdf")
    println("Saved wall comparison plot with different spline orders.")
end