using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
using LsqFit
using CairoMakie: CairoMakie

flow_data = FlowData("examples/data/sod_shock_orb.tape", false)

shock_positions_over_time, angle_estimated = detect(flow_data, ImageProcessingShockDetectionAlgo(0.5, :prewitt))
#create_heatmap_evo_with_shock(flow_data, shock_positions_over_time, angle_estimated, :density_field, true)

dbscan_algo = DBSCANAlgo(0.25, 3, 10)
shock_clusters_over_time = cluster_shock_points(dbscan_algo, shock_positions_over_time, flow_data)
shock_fits_over_time = fit_shock_clusters_over_time(shock_clusters_over_time)

plot_shock_clusters_over_time(shock_clusters_over_time, shock_fits_over_time, flow_data)