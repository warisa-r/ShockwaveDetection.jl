using ShockwaveDetection
include("../NoiseAnalysis.jl")
using .NoiseAnalysis
using Distributions

flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)
noise_data = NoiseAnalysis.NoiseData(0.01, Normal(0, 0.1))
noisy_flow_data = NoiseAnalysis.apply_noise_to_flow(flow_data, noise_data)

point_detect_algo = ImageProcessingShockDetectionAlgo(0.8, :prewitt)
fitting_algo = FittingAlgo(0.1, true)
dbscan_algo = DBSCANAlgo(5.95, 3, 10)

# Detect shock points, clusters, and fits for original and noisy data
original_result = detect(flow_data, point_detect_algo, dbscan_algo, fitting_algo)
noisy_result = detect(noisy_flow_data, point_detect_algo, dbscan_algo, fitting_algo)

rmse_value = NoiseAnalysis.compare_shock_fits_over_time(original_result, noisy_result, flow_data; threshold=10.0)
println("RMSE between fits, cluster_fit_mismatch_count, empty_residuals_count, shock_data_missing_count, residual_length_mismatch_count: ", rmse_value)

#plot_shock_fits_over_time(flow_data, original_detection, true)
#plot_shock_fits_over_time(noisy_flow_data, noisy_detection, true)
