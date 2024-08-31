using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

# Print the current path of this file
println("Current file path: ", @__FILE__)

flow_data = FlowData("data/sod_shock_left_1d.tape", false)

shock_positions_over_time = detect(flow_data, GradientShockDetectionAlgo(0.5))
anim  = create_wave_animation_with_shock(flow_data, shock_positions_over_time)