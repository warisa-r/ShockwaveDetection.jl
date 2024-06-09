# Use Interpolations.jl to perform np.gradient-like operations
using LinearAlgebra: I
using Unitful: ustrip
using Interpolations: interpolate, Gridded, Linear, gradient
using ShockwaveProperties: state_behind, ConservedProps, DRY_AIR, mach_number, shock_temperature_ratio, shock_density_ratio, shock_pressure_ratio

# Function to compute gradients using Interpolations.jl
# without explicit function definition
function compute_gradients(arr::AbstractVector, x)
    # Convert StepRangeLen to LinRange because that's what Interpolations.jl expects
    x_linrange = LinRange(first(x), last(x), length(x))
    itp = interpolate((x_linrange,), arr, Gridded(Linear()))
    grad = only.(gradient.(Ref(itp), x_linrange))
    return grad
end

#TODO: make this work
function rankine_hugoniot(u_values_t, shock_locations)
    # Define normal and tangential vectors for 1D case because state_behind needs it
    n = [1]
    t = [0]

    # flux for the euler equations from ShockwaveProperties.jl's test since it's not exported
    function F(u::ConservedProps)
        v = u.ρv / u.ρ # velocity
        P = pressure(u; gas = DRY_AIR)
        return vcat(u.ρv', (u.ρv .* v' + I * P), (v .* (u.ρE + P))') # stack row vectors
    end

    for shock_location in shock_locations
        u_L = ConservedProps(u_values_t[:, shock_location])
        if mach_number(u_L; gas = DRY_AIR)[1] <= 0.0
            println("Mach number leq 0 at shock location $shock_location and mach number is $(mach_number(u_L; gas = DRY_AIR))")
        else
            println("Mach number gt 0 at shock location $shock_location and mach number is $(mach_number(u_L; gas = DRY_AIR))")
            println(shock_temperature_ratio(mach_number(u_L; gas = DRY_AIR), n; gas = DRY_AIR))
            println(shock_pressure_ratio(mach_number(u_L; gas = DRY_AIR), n; gas = DRY_AIR))
            println(shock_density_ratio(mach_number(u_L; gas = DRY_AIR), n; gas = DRY_AIR))
            u_R = state_behind(u_L, n, t; gas = DRY_AIR)
            # Check Rankine-Hugonoit conditions
            if all(isapprox.(ustrip.(F(u_L) .* n), ustrip.(F(u_R) .* n); atol=1e-6))
                println("Rankine-Hugonoit conditions are satisfied at shock location $shock_location")
            else
                println("Rankine-Hugonoit conditions are not satisfied at shock location $shock_location")
            end
        end
        
    end

end

# Function to detect shock locations at a given time step
function detect_normal_shocks_at_timestep(u_values_t, density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    # Compute first gradients
    density_grad = compute_gradients(density_at_t, x)
    velocity_grad = compute_gradients(velocity_at_t, x)
    pressure_grad = compute_gradients(pressure_at_t, x)

    # Find points where the gradient exceeds a certain threshold
    shock_location_density = findall(gradient -> abs(gradient) > threshold, density_grad)
    shock_location_velocity = findall(gradient -> abs(gradient) > threshold, velocity_grad)
    shock_location_pressure = findall(gradient -> abs(gradient) > threshold, pressure_grad)

    # Combine detections (common shock location across variables)
    shock_locations = intersect(intersect(shock_location_density, shock_location_velocity), shock_location_pressure)

    #Comment this out everything works for now
    #rankine_hugoniot(u_values_t, shock_locations)
    
    return shock_locations
end

"""
    detect_normal_shock(density_field, velocity_field, pressure_field, x0_xmax, t_values; threshold=0.5)

Detects shock positions over time based on given density, velocity, and pressure fields.

# Arguments
- `u_values::Array`: Array containing the conserved state variables at each point x and time t.
- `density_field::Matrix`: Matrix representing the density field over space and time.
- `velocity_field::Matrix`: Matrix representing the velocity field over space and time.
- `pressure_field::Matrix`: Matrix representing the pressure field over space and time.
- `x0_xmax::Tuple`: Tuple containing the start and end positions of the spatial domain.
- `t_values::Vector`: Vector containing time values.
- `threshold::Float64`: Threshold value for detecting shocks. Default is `0.5`.

# Returns
- `shock_positions_over_time::Vector`: Vector containing shock positions at each time step.

# Details
This function iterates over each time step and detects shock positions using the `detect_shocks_at_timestep` function.

"""
function detect_normal_shock(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values; threshold = 0.5)
    len_x = size(density_field, 1)
    x = range(x0_xmax[1], stop=x0_xmax[2], length=len_x)

    shock_positions_over_time = []
    
    for t_step in 1:length(t_values)
        density_field_t = density_field[:, t_step]
        velocity_field_t = velocity_field[:, t_step]
        pressure_field_t = pressure_field[:, t_step]
        u_values_t = u_values[:, :, t_step]
        
        shock_positions = detect_normal_shocks_at_timestep(u_values_t, density_field_t, velocity_field_t, pressure_field_t, x, threshold)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end