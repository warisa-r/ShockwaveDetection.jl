# Use Interpolations.jl to perform np.gradient-like operations
using Interpolations: interpolate, Gridded, Linear, gradient

# Function to compute gradients using Interpolations.jl
# without explicit function definition
function compute_gradients(arr::AbstractVector, x)
    # Convert StepRangeLen to LinRange because that's what Interpolations.jl expects
    x_linrange = LinRange(first(x), last(x), length(x))
    itp = interpolate((x_linrange,), arr, Gridded(Linear()))
    grad = only.(gradient.(Ref(itp), x_linrange))
    return grad
end

# Function to detect shock locations at a given time step
function detect_shocks_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

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
    
    return shock_locations
end

"""
    detect_shock(density_field, velocity_field, pressure_field, x0_xmax, t_values; threshold=0.5)

Detects shock positions over time based on given density, velocity, and pressure fields.

# Arguments
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
function detect_shock(density_field, velocity_field, pressure_field, x0_xmax, t_values; threshold = 0.5)
    len_x = size(density_field, 1)
    x = range(x0_xmax[1], stop=x0_xmax[2], length=len_x)

    shock_positions_over_time = []
    
    for t_step in 1:length(t_values)
        density_field_t = density_field[:, t_step]
        velocity_field_t = velocity_field[:, t_step]
        pressure_field_t = pressure_field[:, t_step]
        
        shock_positions = detect_shocks_at_timestep(density_field_t, velocity_field_t, pressure_field_t, x, threshold)
        push!(shock_positions_over_time, shock_positions)
    end
    
    return shock_positions_over_time
end