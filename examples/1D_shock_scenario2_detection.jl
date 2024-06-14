using ShockwaveDetection
using ShockwaveProperties

flow_data = FlowData("examples/data/euler_scenario_2.out")
detection_alg = EntropyGradientShockDetectionAlgo(DRY_AIR, 0.2)

# Plot the shock create_wave_animation
shock_positions_over_time = detect(flow_data, detection_alg)
anim_with_shock = create_wave_animation_with_shock(flow_data, shock_positions_over_time)