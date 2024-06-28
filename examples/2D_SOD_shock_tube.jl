using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
using LsqFit
using Plots

flow_data = FlowData("examples/data/sod_shock_orb.tape", false)


shock_positions_over_time, angle_estimated = detect(flow_data, ImageProcessingShockDetectionAlgo(0.7, :prewitt))
#create_heatmap_evo_with_shock(flow_data, shock_positions_over_time, angle_estimated, :density_field, true)

dbscan_algo = DBSCANAlgo(0.5, 3, 10)
shock_clusters_over_time = cluster_shock_points(dbscan_algo, shock_positions_over_time, flow_data)
shock_clusters_t = shock_clusters_over_time[250]

shock_cluster = shock_clusters_t[1]

#plot_shock_lines(shock_clusters_t, flow_data)
# Extract x and y coordinates
x_coords = [point[1] for point in shock_cluster]
y_coords = [point[2] for point in shock_cluster]

# Scatter plot for all points in shock_cluster
scatter(x_coords, y_coords, label="Shock Points", marker=:circle)

# Highlight the first point
scatter!([x_coords[1]], [y_coords[1]], label="First Point", color=:red, marker=:star5)

# Highlight the last point
scatter!([x_coords[end]], [y_coords[end]], label="Last Point", color=:green, marker=:hexagon)

# Display the plot
plot!()