using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

# Print the current path of this file
println("Current file path: ", @__FILE__)

# Get the directory of the current script
script_dir = @__DIR__

# Construct the absolute path to the data file
data_file_path = joinpath(script_dir, "data", "sod_shock_left_1d.tape")

# Print the absolute path of the data file
println("Data file path: ", data_file_path)

flow_data = FlowData("data/sod_shock_left_1d.tape", false)

shock_positions_over_time = detect(flow_data, GradientShockDetectionAlgo(0.5))
anim  = create_wave_animation_with_shock(flow_data, shock_positions_over_time)