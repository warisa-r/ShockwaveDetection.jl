using ShockwaveDetection

flow_data = FlowData("examples/data/sod_shock_orb.tape", false)

point_detect_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
dbscan_algo = DBSCANAlgo(0.25, 3, 10)

detection = detect(flow_data, point_detect_algo, dbscan_algo)

plot_shock_clusters_over_time(detection.shock_clusters_over_time, detection.shock_fits_over_time, flow_data)