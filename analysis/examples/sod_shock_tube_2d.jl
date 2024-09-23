using ShockwaveDetection
include("../NoiseAnalysis.jl")
using .NoiseAnalysis
using Distributions
ENV["JULIA_NUM_THREADS"] = "4"

flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)
noise_data = NoiseAnalysis.NoiseData(0.01, Normal(0, 0.1))
noisy_flow_data = NoiseAnalysis.apply_noise_to_flow(flow_data, noise_data)

point_detect_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
fitting_algo = FittingAlgo(0.1, true)
dbscan_algo = DBSCANAlgo(0.25, 3, 10)

original_shock_positions_over_time = ShockwaveDetection.detect_shock_points(flow_data, point_detect_algo)
noisy_shock_positions_over_time = ShockwaveDetection.detect_shock_points(noisy_flow_data, point_detect_algo)

shock_position_difference = NoiseAnalysis.compare_shock_positions_over_time_2d(original_shock_positions_over_time, noisy_shock_positions_over_time)
println("Shock position difference: ", shock_position_difference)
