using ShockwaveDetection

detection_alg = Simple1DShockDetectionAlgo("examples/data/euler_scenario_1.out")

# Plot the animation of the properties over time
#anim = create_wave_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)

# Plot the shock create_wave_animation
shock_positions_over_time = detect(detection_alg, 0.2)
println(shock_positions_over_time)
#anim_with_shock = create_wave_animation_with_shock(x0_xmax, t_values, density_field, velocity_field, pressure_field, shock_positions_over_time)