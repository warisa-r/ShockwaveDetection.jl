using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/supersonic_shock_2.tape")
detection_alg = GradientShockDetectionAlgo(0.2)

# Plot the animation of the properties over time
#anim = create_wave_animation(flow_data)

# Plot the shock create_wave_animation
#shock_positions_over_time = detect(flow_data, detection_alg)
detection = detect(flow_data, detection_alg)
anim_with_shock = create_wave_animation_with_shock(flow_data, detection)
#create_heatmap_evo_with_shock(flow_data, shock_positions_over_time, :density_field)
# Check each member of shock_positions_over_time
