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
function detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo)
    has_obstacle = isnothing(flow_data.u)  # Check if the flow data has an obstacle
    shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo, has_obstacle)
    shock_clusters_over_time = cluster_shock_points(cluster_algo, shock_positions_over_time, flow_data)
    shock_fits_over_time = fit_shock_clusters_over_time(shock_clusters_over_time)
    
    return ShockDetectionResult2D(shock_positions_over_time, shock_clusters_over_time, shock_fits_over_time)
end

# Pipeline for 1D shock detection
function detect(flow_data::FlowData, shock_point_algo::Abstract1DShockDetectionAlgo)
    has_obstacle = isnothing(flow_data.u)  # Check if the flow data has an obstacle
    shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo, has_obstacle)
    
    return ShockDetectionResult1D(shock_positions_over_time)
end
