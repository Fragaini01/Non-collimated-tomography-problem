using Distributed 
folderpath_OWCF = "/home/francesco/Documenti/Master/Blok4/PUK/OWCF/"

cd(folderpath_OWCF)
using Pkg
Pkg.activate(".")
using Dates
using FileIO
using LinearAlgebra
using Equilibrium
using Plots
using Random
using Dierckx
using JLD2

include("Functions.jl")


# Variables
verbose = true
plot = true 

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
radii, zs = regrid_wall(wall , .0001)


if length(ARGS) < 1
    println("Usage: julia Place_detectors.jl <num_detectors>")
    exit(1)
end

num_detectors = parse(Int, ARGS[1])
detectors = detectors_positions(radii, zs, first_wall_size, second_wall_size, num_detectors)
# Add a custom detector at a user-specified position
detectors_on_file(detectors, "detectors(1).csv")
if verbose
    println("Number of detectors: ", length(detectors))
end


if plot
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
    savefig(plt_detectors_rz, "../my_code/output/plots/all_detectors_rz.pdf")
    println("Saved plot of all detectors in R-z plane.")
end