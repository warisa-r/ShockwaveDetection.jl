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
shock_clusters_t = shock_clusters_over_time[70]

shock_cluster = shock_clusters_t[1]

plot_shock_lines(shock_clusters_t, flow_data)