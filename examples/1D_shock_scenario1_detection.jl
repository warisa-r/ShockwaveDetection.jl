using ShockwaveDetection

x0_xmax, t_values, u_values, dims_u = read_output_file("examples/data/euler_scenario_1.out")
density_field, velocity_field, pressure_field = convert_to_primitive(u_values)

# Plot the animation of the properties over time
#anim = create_wave_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)

#Plot the shock create_wave_animation
shock_positions_over_time = detect_normal_shock(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values)
#anim_with_shock = create_wave_animation_with_shock(x0_xmax, t_values, density_field, velocity_field, pressure_field, shock_positions_over_time)