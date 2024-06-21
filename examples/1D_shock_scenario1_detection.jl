using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/supersonic_shock_1.out")
detection_alg = GradientRHShockDetectionAlgo(DRY_AIR, 0.2, 1)

# Plot the animation of the properties over time
#anim = create_wave_animation(flow_data)

# Plot the shock create_wave_animation
shock_positions_over_time = detect(flow_data, detection_alg)
anim_with_shock = create_wave_animation_with_shock(flow_data, shock_positions_over_time)
# Check each member of shock_positions_over_time
