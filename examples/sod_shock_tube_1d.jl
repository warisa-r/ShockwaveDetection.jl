using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
using Distributions

# Create a NoiseData instance
#noise_data = NoiseData(0.01, Normal(0, 1))  # 1% noise intensity, Gaussian distribution

flow_data = FlowData("examples/data/sod_shock_left_1d.tape", false)

point_detect_algo = GradientShockDetectionAlgo(0.5)

detection = detect(flow_data, point_detect_algo)

#plot_shock_fits_over_time(flow_data, detection, true)

shock_positions_over_time = detect(flow_data, point_detect_algo)
anim  = create_wave_animation_with_shock(flow_data, detection)