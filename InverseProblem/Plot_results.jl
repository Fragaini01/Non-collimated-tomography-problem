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
using PyCall
using CSV, DataFrames
using NearestNeighbors
using JLD2
using Measures

include("../Forward_Problem/Functions.jl")

verbose = false

# TRUE VALUES
true_file = "/home/francesco/Documenti/Master/Blok4/PUK/my_code/output/measurements_true_values.jld2"
data_true = load(true_file)
intensities_true = data_true["values"]
D1_true = data_true["D1_array"]
D2_true = data_true["D2_array"]

filepath_equil = "/home/francesco/Documenti/Master/Blok4/PUK/OWCF/equilibrium/ITER/test/80MW_equilibrium.geqdsk"
M, wall = read_geqdsk(filepath_equil, clockwise_phi=false)
radii, zs = regrid_wall(wall , .001)

folder = "/home/francesco/Documenti/Master/Blok4/PUK/my_code/output/inverse/"
n_detectors = 50
files = filter(f -> occursin("$(n_detectors)_W", f) && endswith(f, ".jld2"), readdir(folder; join=true))
println("Trovati $(length(files)) file .jld2 nella cartella $folder")

# Aggiungi il true_file come primo elemento della lista
all_files = vcat([true_file], files)
D1 = D1_true
D2 = D2_true

# Prepara i dati per tutti i file prima di plottare
all_intensities = []
all_titles = []
all_detectors = []
n_plots = length(all_files)
selected_indices = vcat(1, (n_plots-1):3:n_plots)
all_files = [true_file; files[selected_indices[2:end]]]

for (i, file) in enumerate(all_files)
    if file == true_file
        # True values
        intensities = intensities_true
        plot_title = "True Intensities"
    else
        data = load(file)
        # Prendi solo la soluzione migliore (ad esempio la prima)
        new_data = data["sols"]
        n_detectors = size(data["S_abscissa_1"],1)
        push!(all_detectors, n_detectors)
        intensities = new_data[:, :, 12]
        plot_title = "Reconstructed"
    end
    if ndims(intensities) == 1
        intensities = reshape(intensities, length(D1), length(D2))
    end
    push!(all_intensities, intensities)
    push!(all_titles, plot_title)
end

# Calcola il massimo valore assoluto tra tutte le intensità per la colorbar
vmax = maximum([maximum(abs, intensities) for intensities in all_intensities])

psi_mag, psi_bdry = psi_limits(M)
plot_font = "Computer Modern"
Plots.default(fontfamily=plot_font , guidefont = 14, tickfont = 14, legendfont = 14 , titlefont = 20)

@assert length(all_intensities) == length(all_titles) "The number of intensities and titles does not match."
@assert length(all_intensities) == length(all_detectors)+1 "The number of intensities and detectors does not match."
@assert length(all_intensities) == length(selected_indices) "The number of intensities and selected indices does not match."

n_plots = length(all_intensities)
println("Saving data from files")
for i in 1:n_plots
    println(split(all_files[i], "_", limit = 3)[3])
end

plt = plot(layout=(1, length(selected_indices)), aspect_ratio=:equal, size=(600 * length(selected_indices), 850))
color = cgrad([:white, :blue, :red], [0.0, 0.5, 1.0])


clims_common = (0.0, maximum(intensities_true))
xlims = (min(minimum(D1), minimum(wall.r)), max(maximum(D1), maximum(wall.r)))
ylims = (min(minimum(D2), minimum(wall.z)), max(maximum(D2), maximum(wall.z)))

ticks = [0, 2.5e20, 5e20, 7.5e20, 1e21]
ticklabels = ["0", "2.5·10²⁰", "5·10²⁰", "7.5·10²⁰", "10²¹"]

for i in 1:n_plots
    intensities = all_intensities[i]
    title = all_titles[i]
    scatter!(
        plt[i],
        repeat(D1, inner=length(D2)),
        repeat(D2, outer=length(D1)),
        marker_z = vec(intensities),
        color = color,
        clims = clims_common,
        ms = 6,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        markershape = :circle,
        title = title,
        xlabel = "R [m]",
        ylabel = "Z [m]",
        colorbar = (i == 2),
        colorbar_title = "I [#neutrons/voxel volume/s]",
        colorbar_titlefont = 14,
        colorbar_tickfont = 14,
        colorbar_orientation = :vertical,
        colorbar_position = (i == 2 ? :left : :none),
        colorbar_ticks = (ticks, ticklabels),
        aspect_ratio = :equal,
        label = "",
        xlims = xlims,
        ylims = ylims
    )
    plot!(plt[i], radii, zs, label="", linecolor=:black, linewidth=1)
    scatter!(plt[i], wall.r, wall.z, label="", linewidth=2.5, color=:black)
end

# Salva il plot unico
savefig(plt, joinpath(folder, "plot_all_results_$(n_detectors).pdf"))


last_file = files[end]

data = load(last_file)


new_data = data["sols"]
n_detectors = size(data["S_abscissa_1"],1)


solution_numbers = 2:4:20
n_solutions = length(solution_numbers)
n_rows = 2
n_cols = ceil(Int, n_solutions / n_rows)
plt_crs = plot(
    layout = (n_rows, n_cols),
    #aspect_ratio = :equal,
    size = (400 * n_cols, 450 * n_rows),  
    background_color = :white,
    legend = false,
    grid = false,
    framestyle = :box,
    left_margin = 10mm,
    right_margin = 10mm,
    top_margin = 10mm,
    bottom_margin = 10mm,
    subplot_padding = 10mm
)

l_curve_x = data["l_curve_x"]
l_curve_y = data["l_curve_y"]
#clims_common = (0.0, maximum(intensities_true))
for (i, solution_number) in enumerate(solution_numbers)
    intensities = new_data[:, :, solution_number]
    if ndims(intensities) == 1
        intensities = reshape(intensities, length(D1), length(D2))
    end

    scatter!(
        plt_crs[i],
        repeat(D1, inner=length(D2)),
        repeat(D2, outer=length(D1)),
        marker_z = vec(intensities),
        color = color,
        clims = clims_common,
        ms = 6,
        markerstrokecolor = :transparent,
        markerstrokewidth = 0,
        markershape = :circle,

        colorbar = false,
        colorbar_title = "I [#neutrons/voxel volume/s]",
        colorbar_titlefont = 12,
        colorbar_tickfont = 12,
        colorbar_orientation = :vertical,
        colorbar_position = :right,
        
        label = "",
        xlims = xlims,  
        ylims = ylims,   
        aspect_ratio = :equal,
        xlabel = "R [m]",
        ylabel = "Z [m]",
    )
    plot!(plt_crs[i], radii, zs, label="", linecolor=:black, linewidth=1)
    scatter!(plt_crs[i], wall.r, wall.z, label="", linewidth=2.5, color=:black)
    plot!(plt_crs[i], xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16)
    plot!(plt_crs[i], legend=false)
    plot!(plt_crs[i], title="λ = $(round(l_curve_y[solution_number]; digits=2))", titlefontsize=14)

    if i == n_solutions
        plot!(
            plt_crs[i+1],
            l_curve_x,
            l_curve_y,
            linecolor = :orange,
            linestyle = :solid,
            #marker = :circle,
            linewidth = 2,
            label = "L-curve",
            xlabel = "log(||Ax-b||)",
            ylabel = "log(||Lx||)",
            xguidefontsize = 16,
            yguidefontsize = 16,
            legend = false,
            xlims = (minimum(l_curve_x) - 0.5, maximum(l_curve_x) + 0.5),
            ylims = (-0.01, maximum(l_curve_y) + 0.01),
            title = "L-curve",
            titlefontsize = 14,
            xtickfontsize = 14,
            ytickfontsize = 14
        )
    end
end

savefig(plt_crs, joinpath(folder, "plot_solutions_subplot_$(n_detectors).pdf"))
