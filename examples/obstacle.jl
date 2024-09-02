using ShockwaveDetection
using ShockwaveProperties
using Euler2D:Euler2D

flow_data = FlowData("examples/data/obstacle/circular_obstacle_radius_1.celltape")
point_detect_algo = ImageProcessingShockDetectionAlgo(0.2, :prewitt)
dbscan_algo = DBSCANAlgo(0.25, 3, 10)

detection = detect(flow_data, point_detect_algo, dbscan_algo)

#plot_shock_fits_over_time(flow_data, detection, false)
create_heatmap_evo_with_shock(flow_data, detection, :density_field, true, false)
