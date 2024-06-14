using ShockwaveProperties

# Upstream (Approaching the Shock): The entropy gradient typically decreases or remains constant.
# Across the Shock: There is a sudden increase in entropy, causing the entropy gradient to change abruptly.
# Downstream (After the Shock): The entropy gradient stabilizes at a higher value due to the increased entropy.
# We assume that at the boundary, there is no shock (kind of true looking at the animated plots) 
function shock_entropy_change(discontinuity_locations, density_at_t, pressure_at_t, x, gas)
    c_v = gas.c_v
    gamma = gas.γ

    # TODO: check proper source, not just use the pic u see in reddit https://www.reddit.com/r/thermodynamics/comments/113vvxq/assuming_ideal_gas_and_caloricallyperfect/
    # Ignore the constant because it's going to get canceled out in the gradient calculation
    entropy_at_t = c_v * log.(pressure_at_t ./ density_at_t.^gamma)

    entropy_grad = compute_gradients(entropy_at_t, x)

    # Find locations where entropy gradient changes sign (indicating shock)
    shock_location_entropy = findall(gradient -> sign(gradient) != sign(entropy_grad[1]), entropy_grad)

    shock_locations = intersect(discontinuity_locations, shock_location_entropy)

    return shock_locations
end

struct EntropyGradientShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    gas::CaloricallyPerfectGas
    threshold::Float64
end # EntropyGradientShockDetectionAlgo

function detect_entropy_normal_shocks_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, alg)
    threshold = alg.threshold
    gas = alg.gas

    # Detect discontinuities in the flow at a given time step
    discontinuity_locations = detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    #  Check if Rankine-Hugoniot conditions are satisfied
    shock_locations = shock_entropy_change(discontinuity_locations, density_at_t, pressure_at_t, x, gas)
    
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
        shock_positions = detect_entropy_normal_shocks_at_timestep(density_field_t, velocity_field_t, pressure_field_t, x, alg)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end
