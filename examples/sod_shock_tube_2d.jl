using ShockwaveDetection
ENV["JULIA_NUM_THREADS"] = "4"

flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)

point_detect_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
fitting_algo = FittingAlgo(0.1, true)
dbscan_algo = DBSCANAlgo(0.25, 3, 10)

detection = detect(flow_data, point_detect_algo, dbscan_algo, fitting_algo)

plot_shock_fits_over_time(flow_data, detection, true)
#create_heatmap_evo_with_shock(flow_data, detection, :density_field, true, false)