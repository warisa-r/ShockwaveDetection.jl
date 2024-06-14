# Use Interpolations.jl to perform np.gradient-like operations
using LinearAlgebra: I
using Unitful: ustrip
using ShockwaveProperties: state_behind, ConservedProps, DRY_AIR, mach_number, shock_temperature_ratio, shock_density_ratio, shock_pressure_ratio



abstract type AbstractShockDetectionAlgo end
# Abstract base type for 1D shock detectors
abstract type Abstract1DShockDetectionAlgo <: AbstractShockDetectionAlgo end

#TODO: more efficient?
function shock_rankine_hugoniot(discontinuity_locations, u_values_t)
    # Define normal and tangential vectors for 1D case because state_behind needs it
    n = [1]
    t = [0]

    normal_shock_positions = []

    # Flux for the euler equations from ShockwaveProperties.jl's test since it's not exported
    function F(u::ConservedProps)
        v = u.ρv / u.ρ # velocity
        P = pressure(u; gas = DRY_AIR)
        return vcat(u.ρv', (u.ρv .* v' + I * P), (v .* (u.ρE + P))') # stack row vectors
    end

    for discontinuity in discontinuity_locations
        
        u_front = ConservedProps(u_values_t[:, discontinuity])
        # Check if Mach number is greater than 1 at the shock location -> for normal shock in supersonic flow
        # else, the shock is not a normal shock and this is not considered
        if abs(mach_number(u_front; gas = DRY_AIR)[1]) >= 1
            # Function implemented in ShockwaveProperties.jl from Anderson&Anderson(only holds for Mach number >= 1)
            # But in reality the Rankine-Hugonoit conditions hold for all discontinuities
            u_behind = state_behind(u_front, n, t; gas = DRY_AIR)
            # Check Rankine-Hugonoit conditions
            if all(isapprox.(ustrip.(F(u_front) .* n), ustrip.(F(u_behind) .* n); atol=1e-6))
                println("Rankine-Hugonoit conditions are satisfied at shock location $discontinuity")
                push!(normal_shock_positions, discontinuity)
            else
                # Debug message
                println("Rankine-Hugonoit conditions are not satisfied at shock location $discontinuity, mach number = $(mach_number(u_front; gas = DRY_AIR))")
            end
        else
            println("Mach number is less than 1")
        end
    end
    
    return normal_shock_positions
end

# TODO: Add gas to this struct
# RHS stands for Rankine-Hugoniot Shock
struct GradientRHShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    threshold::Float64
end # GradientRHShockDetectionAlgo

# This function detects shock locations in a fluid at a given time step. It computes the gradients of density, velocity, and pressure, 
# and identifies locations where these gradients exceed a threshold (scaled by the maximum gradient of each variable). 
# It then checks if the Rankine-Hugoniot conditions are satisfied at these locations.
function detect_RHS_normal_shocks_at_timestep(u_values_t, density_at_t, velocity_at_t, pressure_at_t, x, alg)
    threshold = alg.threshold

    # Detect discontinuities in the flow at a given time step
    discontinuity_locations = detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    #  Check if Rankine-Hugoniot conditions are satisfied
    shock_locations = shock_rankine_hugoniot(discontinuity_locations, u_values_t)
    
    return shock_locations
end

function detect(flow_data::FlowData, alg::GradientRHShockDetectionAlgo)
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
        shock_positions = detect_RHS_normal_shocks_at_timestep(u_values_t, density_field_t, velocity_field_t, pressure_field_t, x, alg)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end
