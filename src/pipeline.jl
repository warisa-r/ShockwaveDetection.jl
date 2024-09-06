abstract type AbstractShockDetectionResult end

struct ShockDetectionResult2D <: AbstractShockDetectionResult
    shock_positions_over_time::Vector{Any}
    shock_clusters_over_time::Vector{Any}
    shock_fits_over_time::Vector{Any}
end

struct ShockDetectionResult1D <: AbstractShockDetectionResult
    shock_positions_over_time::Vector{Any}
end

# Pipeline for 2D shock detection
"""
    detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo)

Detects shocks in 2D flow data, clusters the shockpoints and applies fitting to the cluster.

# Arguments
- `flow_data::FlowData`: A `FlowData` object containing the 2D flow field data.
- `shock_point_algo::Abstract2DShockDetectionAlgo`: An algorithm for detecting shock points in 2D flow data.
- `cluster_algo::DBSCANAlgo`: A clustering algorithm (e.g., DBSCAN) to group detected shock points into clusters.

# Returns
- `ShockDetectionResult2D`: An object containing:
  - `shock_positions_over_time`: Detected shock points over time.
  - `shock_clusters_over_time`: Clusters of shock points over time.
  - `shock_fits_over_time`: Fitted shock curves over time.

# Description
This function detects shock points in 2D flow data using a specified shock detection algorithm. Detected shock points are clustered using the provided `DBSCANAlgo`, and then the clusters are fitted to create a smooth representation of the shock over time.
"""
function detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo)
    has_obstacle = isnothing(flow_data.u)  # Check if the flow data has an obstacle
    shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo, has_obstacle)
    shock_clusters_over_time = cluster_shock_points(cluster_algo, shock_positions_over_time, flow_data)
    shock_fits_over_time = fit_shock_clusters_over_time(shock_clusters_over_time)
    
    return ShockDetectionResult2D(shock_positions_over_time, shock_clusters_over_time, shock_fits_over_time)
end

"""
    detect(flow_data::FlowData, shock_point_algo::Abstract1DShockDetectionAlgo)

Detects shocks in 1D flow data.

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
    has_obstacle = isnothing(flow_data.u)  # Check if the flow data has an obstacle
    shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo, has_obstacle)
    
    return ShockDetectionResult1D(shock_positions_over_time)
end
