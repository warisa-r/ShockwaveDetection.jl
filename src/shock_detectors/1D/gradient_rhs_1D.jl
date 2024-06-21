using Unitful: ustrip
using ShockwaveProperties: ShockwaveProperties

#TODO: more efficient?
function shock_rankine_hugoniot(discontinuity_locations, u_values_t, gas, hugoniot_tol=1e-6)
    # Define normal and tangential vectors for 1D case because state_behind needs it
    n = [1]
    t = [0]

    normal_shock_positions = []

    for discontinuity in discontinuity_locations
        # Check if discontinuity is at the first or last point of the grid
        if discontinuity < 10 || discontinuity >= size(u_values_t, 2)-10
            println("Discontinuity is at the edge of the grid at index $discontinuity.")
            continue
        end
        # Look for the state behind and in front of the shock
        # 10 grid points are how far we are looking (hard-coded)
        u1 = ShockwaveProperties.ConservedProps(u_values_t[:, discontinuity-10])
        u2 = ShockwaveProperties.ConservedProps(u_values_t[:, discontinuity+10])

        # Don't include kinetic energy since the conditions separately account for the kinetic energy through the terms involving velocity
        e1 = ShockwaveProperties.specific_static_internal_energy(u1; gas)
        e2 = ShockwaveProperties.specific_static_internal_energy(u2; gas)

        p1 = ShockwaveProperties.pressure(u1; gas)
        p2 = ShockwaveProperties.pressure(u2; gas)

        ρ1 = ShockwaveProperties.density(u1)
        ρ2 = ShockwaveProperties.density(u2)

        # Hugoniot equation (7.9) from Anderson&Anderson
        if isapprox(ustrip(e2 - e1), ustrip((p1 + p2) * (1/ρ1 - 1/ρ2) / 2), atol=hugoniot_tol)
            push!(normal_shock_positions, discontinuity)
        else
            println("internal energy difference: ", e2 - e1)
            println("righ hand side: ", (p1 + p2) * (1/ρ1 - 1/ρ2) / 2)
        end
    end
    
    return normal_shock_positions
end

# RH stands for Rankine-Hugoniot Shock
struct GradientRHShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    gas::CaloricallyPerfectGas
    threshold::Float64
    hugoniot_tol::Float64
end # GradientRHShockDetectionAlgo

# This function detects shock locations in a fluid at a given time step. It computes the gradients of density, velocity, and pressure, 
# and identifies locations where these gradients exceed a threshold (scaled by the maximum gradient of each variable). 
# It then checks if the Rankine-Hugoniot conditions are satisfied at these locations.
function detect_RHS_normal_shocks_at_timestep(u_values_t, density_at_t, velocity_at_t, pressure_at_t, x, alg)
    threshold = alg.threshold
    gas = alg.gas
    hugoniot_tol = alg.hugoniot_tol

    # Detect discontinuities in the flow at a given time step
    discontinuity_locations = detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    #  Check if Rankine-Hugoniot conditions are satisfied
    shock_locations = shock_rankine_hugoniot(discontinuity_locations, u_values_t, gas, hugoniot_tol)
    
    return shock_locations
end

function detect(flow_data::FlowData, alg::GradientRHShockDetectionAlgo)
    # Unpack all the values from the detector
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field
    bounds = flow_data.bounds
    tsteps = flow_data.tsteps
    u = flow_data.u

    len_x = size(density_field, 1)
    x = range(bounds[1][1], stop=bounds[1][2], length=len_x)

    shock_positions_over_time = []
    
    for t_step in 1:length(tsteps)
        density_field_t = density_field[:, t_step]
        velocity_field_t = velocity_field[:, t_step]
        pressure_field_t = pressure_field[:, t_step]
        u_values_t = u[:, :, t_step]
        
        # Use the simple shock detection algorithm to detect the normal shock
        shock_positions = detect_RHS_normal_shocks_at_timestep(u_values_t, density_field_t, velocity_field_t, pressure_field_t, x, alg)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end
