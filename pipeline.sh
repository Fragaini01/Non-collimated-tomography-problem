#!/bin/bash
rm output/*.vc
rm output/detector_*.png
rm output/a_detector_*.png
julia Forward_Problem/Place_detectors.jl 99
julia Forward_Problem/Multi_detector_grids_parallel.jl
julia Forward_Problem/compute_a.jl
julia Forward_Problem/compute_Intensities.jl