# ShockDetection.jl

ShockDetection.jl is a Julia package for detecting shocks in time series data from Godunov solver.
The packages can detect shock discontinuities for grid data with dimensions
- Nu x Nx x Nt when Nx are the number of grid points in x axis and Nt is the number of time steps in the simulation
- Nu x Nx x Ny x Nt when Nx and Ny are the number of grid points in x and y axes respectively, and Nt is the number of time steps in the simulation

The package can only process simulation files that end with .tape and .celltape. Other simulation file types are not supported!

## Features
- **Shock detection (1D and 2D)**: 
  - For 1D data (Nu x Nx x Nt): Detects shocks based on gradients using the finite difference method. 
  - For 2D data (Nu x Nx x Ny x Nt): Detects shocks based on gradients, clusters points using DBSCAN and fits optimal shock curves via the Levenberg–Marquardt algorithm.
- **Visualization**: Generate shock visualizations, calculate normal vectors for 2D simulations with `create_heatmap_evo_with_shock`, or use `calculate_normal_vector` to compute vectors at specific time frames.
- **Compatibility**: Works with Julia 1.9 and 1.10.

## Installation

You can install ShockDetection.jl using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode
There are several dependencies needed to run this pacakge. Check out `Project.toml` to see what packages you need to have in order to run `ShockwaveDetection.jl`

## Guide

Go to our examples directory to see how to use our functions!

### Short guide on how to use our package on grid data with dimension Nx x Nt
When you wish to use our package to detect shock discontinuities,  you should first input your data file path into `FlowData`.
After that, define your algorithm of finite difference that best suited your simulation data and use the function `detect` to automate the detection process. The results of your detection can be found in `ShockDetectionResult1D`. You can further use various functions in visualization to see the result.

    Example:
    ```julia
    flow_data = FlowData("examples/data/supersonic_shock_2.tape")
    detection_algo = GradientShockDetectionAlgo(threshold=0.5)
    result = detect(flow_data, detection_algo)

### Short guide on how to use our package on grid data with dimension Nx x Ny x Nt
When you wish to use our package to detect shock discontinuities,  you should first input your data file path into `FlowData`.
After that, define your algorithm of finite difference, clustering, and fitting that best suited your simulation data and use
the function `detect` to automate the detection process. The results of your detection can be found in `ShockDetectionResult2D`. You can further use various functions in visualization to
see the result.

    Example:
    ```julia
    flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)
    detection_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
    fitting_algo = FittingAlgo(0.1, true)
    dbscan_algo = DBSCANAlgo(0.25, 3, 10)
    result = detect(flow_data, detection_algo, dbscan_algo, fitting_algo)

# Functions and structs

## Main functions and structs of our package
These are the most important functions and structs in our package. Others are functions/structs called/used/stored by `detect`/ `ShockDetectionResult1D` and `ShockDetectionResult2D` to perform shock detections but they still exist in this documentation to clarify to you how our algorithm works!
```@docs
detect
ShockDetectionResult1D
ShockDetectionResult2D
```

## Input-related structs and functions
```@docs
convert_to_primitive
FlowData
```

## Finite difference-related functions and structs
```@docs
ImageProcessingShockDetectionAlgo
GradientShockDetectionAlgo
```

## Cluster-related functions and structs
```@docs
cartesian_index_to_xy
cluster_shock_points
DBSCANAlgo
```

## Fitting-related functions and structs
```@docs
FittingAlgo
fit_shock_clusters_over_time
calculate_normal_vector
```

## Visualization functions
```@docs
create_heatmap_evo_with_shock
create_wave_animation
create_wave_animation_with_shock
plot_shock_fits_over_time
```