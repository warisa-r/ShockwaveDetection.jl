using ShockwaveProperties

function shock_entropy_change(discontinuity_locations, u_values_t, gas = DRY_AIR)
    shock_locations = []
    
    return shock_locations
end

struct EntropyGradientShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    threshold::Float64
    gas::CaloricallyPerfectGas
end # EntropyGradientShockDetectionAlgo

function detect_entropy_normal_shocks_at_timestep(u_values_t, density_at_t, velocity_at_t, pressure_at_t, x, alg)
    threshold = alg.threshold
    gas = alg.gas

    # Detect discontinuities in the flow at a given time step
    discontinuity_locations = detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    #  Check if Rankine-Hugoniot conditions are satisfied
    shock_locations = shock_entropy_change(discontinuity_locations, u_values_t, gas)
    
    return shock_locations
end

function detect(flow_data::FlowData, alg::EntropyGradientShockDetectionAlgo)
    # Unpack all the values from the detector
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field
    x0_xmax = flow_data.x0_xmax
    t_values = flow_data.t_values
    u_values = flow_data.u_values

    len_x = size(density_field, 1)
    x = range(x0_xmax[1], stop=x0_xmax[2], length=len_x)

    shock_positions_over_time = []
    
    for t_step in 1:length(t_values)
        density_field_t = density_field[:, t_step]
        velocity_field_t = velocity_field[:, t_step]
        pressure_field_t = pressure_field[:, t_step]
        u_values_t = u_values[:, :, t_step]
        
        # Use the simple shock detection algorithm to detect the normal shock
        shock_positions = detect_entropy_normal_shocks_at_timestep(u_values_t, density_field_t, velocity_field_t, pressure_field_t, x, alg)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end
