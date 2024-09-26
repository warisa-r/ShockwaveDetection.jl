using ShockwaveDetection
include("../NoiseAnalysis.jl")
using .NoiseAnalysis
using Distributions
ENV["JULIA_NUM_THREADS"] = "4"

flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)
noise_data = NoiseAnalysis.NoiseData(0.01, Normal(0, 0.1))
noisy_flow_data = NoiseAnalysis.apply_noise_to_flow(flow_data, noise_data)

point_detect_algo = ImageProcessingShockDetectionAlgo(0.8, :prewitt)
fitting_algo = FittingAlgo(0.1, true)
dbscan_algo = DBSCANAlgo(1.05, 3, 10)

#original_detection = detect(flow_data, point_detect_algo, dbscan_algo, fitting_algo)
#noisy_detection = detect(noisy_flow_data, point_detect_algo, dbscan_algo, fitting_algo)

original_shock_positions_over_time = ShockwaveDetection.detect_shock_points(flow_data, point_detect_algo)
noisy_shock_positions_over_time = ShockwaveDetection.detect_shock_points(noisy_flow_data, point_detect_algo)
#println("original_shock_positions_over_time: ", original_shock_positions_over_time)
#println("noisy_shock_positions_over_time: ", noisy_shock_positions_over_time)

original_shock_clusters_over_time = cluster_shock_points(dbscan_algo, original_shock_positions_over_time, flow_data)
noisy_shock_clusters_over_time = cluster_shock_points(dbscan_algo, noisy_shock_positions_over_time, flow_data)
#println("original_shock_clusters_over_time: ", original_shock_clusters_over_time)
#println("noisy_shock_clusters_over_time: ", noisy_shock_clusters_over_time)

original_shock_fits_over_time = fit_shock_clusters_over_time(original_shock_clusters_over_time, fitting_algo)
noisy_shock_fits_over_time = fit_shock_clusters_over_time(noisy_shock_clusters_over_time, fitting_algo)
#println("original_shock_fits_over_time: ", original_shock_fits_over_time)
#println("noisy_shock_fits_over_time: ", noisy_shock_fits_over_time)


rmse_value = NoiseAnalysis.compare_shock_fits_over_time(original_shock_fits_over_time, noisy_shock_fits_over_time)
println("RMSE between fits: ", rmse_value)

#shock_fits_difference = NoiseAnalysis.compare_shock_fits_over_time(original_shock_fits_over_time, noisy_shock_fits_over_time)

#plot_shock_fits_over_time(flow_data, original_detection, true)
#plot_shock_fits_over_time(noisy_flow_data, noisy_detection, true)
#shock_position_difference = NoiseAnalysis.compare_shock_positions_over_time_2d(original_shock_positions_over_time, noisy_shock_positions_over_time)
#println("Shock position difference: ", shock_position_difference)
