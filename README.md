# Non-Collimated Tomography Problem

Repository for the PUK project at the University of Copenhagen

This repository contains the core code to perform both the forward model and the inverse problem for a set of detectors in a non-collimated configuration, as installed on ITER.

Note: All code depends on the `OWCF` folder, which is required for execution.

-------------------------------------------------------------------------------

## Repository Structure

### pipeline.sh
A bash script that runs the entire forward model pipeline.

- Default number of detectors: 99 (can be modified as needed).

-------------------------------------------------------------------------------

## Forward Problem

Located in the `Forward_Problem/` folder.

### Place_detectors.jl
Places a specified number of detectors and outputs a file `detectors.csv` containing their coordinates.

Usage:
    julia Place_detectors.jl 

### Multi_detectors_grids_parallel.jl
Computes the voxels seen by each detector using the `createCustomLOS.jl` script from the `OWCF` folder. Generates a `.vc` file for each detector. The number of workers for parallel computation can be set within the script.

### compute_a.jl
Generates matrix A, which encodes the physical constraints of the forward problem. Outputs the matrix in a file named `a.jl`.

### compute_intensities.jl
Simulates intensity measurements and produces the following output files:

- measurements_true_values.jld2: true voxel intensities.
- measurements_20.jld2, measurements_50.jld2, measurements.jld2: simulated measurements for 20, 50, and 99 detectors, respectively.

### plot_for_report.jl
Generates four plots to visualize the forward model:

1. Voxels seen by each detector in the x-y plane.
2. Each row of matrix A.
3. Comparison of spline interpolations for the wall.
4. Sum of all rows in matrix A.

Requires the relevant input files to be provided.

-------------------------------------------------------------------------------

## Inverse Problem

Located in the `Inverse_Problem/` folder.

### start_solveInverseProblem.jl
Launches the inverse problem solver, based on the `OWCF` framework.

### Plot_results.jl
Allows selection of the number of detectors and generates two plots:

1. Comparison between true intensities and the best reconstruction.
2. Comparison of four reconstructions using different lambda values.
