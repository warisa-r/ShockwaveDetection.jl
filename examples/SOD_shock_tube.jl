using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/sod1.out", false)

create_tube_field_evo(flow_data, :velocity_field)

detection_alg = GradientEntropyShockDetectionAlgo(DRY_AIR, 0.5)

# Plot the shock create_wave_animation
shock_positions_over_time = detect(flow_data, detection_alg)
#anim_with_shock = create_wave_animation_with_shock(flow_data, shock_positions_over_time)

create_tube_field_evo_with_shock(flow_data, shock_positions_over_time, :velocity_field)

