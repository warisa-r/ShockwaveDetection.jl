using TimerOutputs

abstract type AbstractShockDetectionResult end

"""
    struct ShockDetectionResult2D <: AbstractShockDetectionResult

A structure to hold the results of 2D shock detection.

# Fields
- `shock_positions_over_time::Vector{Any}`: A vector of the CartesianIndices of shock points detected in every time frame.
- `shock_clusters_over_time::Vector{Any}`: A vector of groups of clustered shock points detected in every time frame.
- `shock_fits_over_time::Vector{Any}`: Fitted shock curves over time. In each time frame is a vector of objects `Fitting`

# Description
This structure is used to store the results of shock detection in 2D flow data. It contains the detected shock points, clusters of shock points, and fitted shock curves over time.
"""
struct ShockDetectionResult2D <: AbstractShockDetectionResult
    shock_positions_over_time::Vector{Any}
    shock_clusters_over_time::Vector{Any}
    shock_fits_over_time::Vector{Any}
end

"""
    struct ShockDetectionResult1D <: AbstractShockDetectionResult

A structure to hold the results of 1D shock detection.

# Fields
- `shock_positions_over_time::Vector{Any}`: A vector of the CartesianIndices of shock points detected in every time frame.

# Description
This structure is used to store the results of shock detection in 1D flow data. It contains the detected shock points over time.
"""
struct ShockDetectionResult1D <: AbstractShockDetectionResult
    shock_positions_over_time::Vector{Any}
end

"""
    detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo, fitting_algo::FittingAlgo)

Detects shocks in 2D flow data, clusters the shockpoints and applies fitting to the cluster. It also shows the runtime and memory allocations required in each subprocess

# Arguments
- `flow_data::FlowData`: A `FlowData` object containing the 2D flow field data.
- `shock_point_algo::Abstract2DShockDetectionAlgo`: An algorithm for detecting shock points in 2D flow data.
- `cluster_algo::DBSCANAlgo`: A clustering algorithm (e.g., DBSCAN) to group detected shock points into clusters.
- `fitting_algo::FittingAlgo`: An algorithm for fitting shock clusters to create a smooth representation of the shock over time.

# Returns
- `ShockDetectionResult2D`: An object containing:
  - `shock_positions_over_time`: Detected shock points over time.
  - `shock_clusters_over_time`: Clusters of shock points over time.
  - `shock_fits_over_time`: Fitted shock curves over time.

# Description
This function detects shock points in 2D flow data using a specified shock detection algorithm. Detected shock points are clustered using the provided `DBSCANAlgo`, and then the clusters are fitted to create a smooth representation of the shock over time.
"""
function detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo, fitting_algo::FittingAlgo)
    to = TimerOutput()
    
    @timeit to "Detect Shock Points(2D)" begin
        shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo)
    end
    
    @timeit to "Cluster Shock Points" begin
        shock_clusters_over_time = cluster_shock_points(cluster_algo, shock_positions_over_time, flow_data)
    end
    
    @timeit to "Fit Shock Clusters" begin
        shock_fits_over_time = fit_shock_clusters_over_time(shock_clusters_over_time, fitting_algo)
    end
    
    show(to, sortby = :firstexec)
    
    return ShockDetectionResult2D(shock_positions_over_time, shock_clusters_over_time, shock_fits_over_time)
end

"""
    detect(flow_data::FlowData, shock_point_algo::Abstract1DShockDetectionAlgo)

Detects shocks in 1D flow data and show the runtime and memory allocations required in each subprocess.

# Arguments
- `flow_data::FlowData`: A `FlowData` object containing the 1D flow field data.
- `shock_point_algo::Abstract1DShockDetectionAlgo`: An algorithm for detecting shock points in 1D flow data.

# Returns
- `ShockDetectionResult1D`: An object containing:
  - `shock_positions_over_time`: Detected shock points over time.

# Description
This function detects shock points in 1D flow data using a specified shock detection algorithm. 
"""
function detect(flow_data::FlowData, shock_point_algo::Abstract1DShockDetectionAlgo)
    to = TimerOutput()

    has_obstacle = isnothing(flow_data.u)  # Check if the flow data has an obstacle

    @timeit to "Detect Shock Points (1D)" begin
        shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo)
    end

    show(to, sortby = :firstexec)
    
    return ShockDetectionResult1D(shock_positions_over_time)
end
