#TODO: pipeline 1D shock detection too

struct ShockDetectionResult2D
    shock_positions_over_time::Vector{Any}
    shock_clusters_over_time::Vector{Any}
    shock_fits_over_time::Vector{Any}
end

# Pipeline is only for 2D shock detection now
function detect(flow_data::FlowData, shock_point_algo::Abstract2DShockDetectionAlgo, cluster_algo::DBSCANAlgo, make_initial_guess::Bool)
    has_obsticle = false
    if isnothing(flow_data.u)
        has_obsticle = true
    end
    
    shock_positions_over_time = detect_shock_points(flow_data, shock_point_algo, has_obstacle)
    shock_clusters_over_time = cluster_shock_points(cluster_algo, shock_positions_over_time, flow_data)
    shock_fits_over_time = fit_shock_clusters_over_time(shock_clusters_over_time, make_initial_guess)
    return ShockDetectionResult2D(shock_positions_over_time, shock_clusters_over_time, shock_fits_over_time)
end
