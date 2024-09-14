# ShockDetection.jl

ShockDetection.jl is a Julia package for detecting shocks in time series data from Godunov solver.
The packages can detect shock discontinuities for grid data with dimensions
- Nx x Nt when Nx are the number of grid points in x axis and Nt is the number of time steps in the simulation
- Nx x Ny x Nt when Nx and Ny are the number of grid points in x and y axes respectively, and Nt is the number of time steps in the simulation

The package can only process simulation files that end with .tape and .celltape. Other simulation file types are not supported!

## Features

- A method for shock detection. The gradient of each point in the simulation is calculated using finite difference method, clustered by DBSCAN algorithm and optimal shock curves are fitted using Levenberg–Marquardt algorithm.
- Visualization is also available within the package. You can also automatically calculate the normal vectors of 2D simulation as you call the function `create_heatmap_evo_with_shock`. But if you wish to obtain the normal vector of each curve in a specified time frame, the function `calculate_normal_vector` is also available for you.
- Compatible with Julia 1.9 and 1.10.

## Installation

You can install ShockDetection.jl using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode
There are several dependencies needed to run this pacakge. Checkout `Project.toml` to see what packages you need to have in order to run `ShockwaveDetection.jl`

## Guide

Go to our examples directory to see how to use our functions!

### Short guide on how to use our package on grid data with dimension Nx x Ny x Nt
When you wish to use our package to detect shock discontinuities,  you should first input your data file path into `FlowData`.
After that, define your algorithm of finite difference, clustering, and fitting that best suited your simulation data and use
the function `detect` to automate the process the detection process. The results of your detection can be found in `ShockDetectionResult2D`. You can further use various functions in visualization to
see the result.

# Functions and structs

## Main functions and structs of our package
These are the most important functions and structs in our package. Others are functions/structs called/used/stored by `detect`/ `ShockDetectionResult2D` to perform shock detections but they still exist in this documentation to clarify to you how our algorithm works!
```@docs
detect
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
```

## Cluster-related functions and structs
```@docs
cartesian_index_to_xy
cluster_shock_points
DBSCANAlgo
```

## Fitting-related functions and structs
```@docs
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