using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
include("../analysis/NoiseAnalysis.jl")
using .NoiseAnalysis
using Distributions


# Create a NoiseData instance
noise_data = NoiseAnalysis.NoiseData(0.01, Normal(0, 1))  # 1% noise intensity, Gaussian distribution

flow_data = FlowData("examples/data/sod_shock_left_1d.tape", false)

# Apply noise to flow data
noisy_flow_data = NoiseAnalysis.apply_noise_to_flow(flow_data, noise_data)

point_detect_algo = GradientShockDetectionAlgo(0.5)

#detection = detect(flow_data, point_detect_algo)

#plot_shock_fits_over_time(flow_data, detection, true)

#shock_positions_over_time = detect(flow_data, point_detect_algo)
original_shock_positions = detect(flow_data, point_detect_algo)
noisy_shock_positions = detect(noisy_flow_data, point_detect_algo)
#println("original_shock_positions: ", original_shock_positions)
#println("noisy_shock_positions: ", noisy_shock_positions)
#anim  = create_wave_animation_with_shock(flow_data, detection)

# Detect and compare shock positions with and without noise
shock_diff = NoiseAnalysis.compare_shock_positions_over_time(original_shock_positions, noisy_shock_positions)