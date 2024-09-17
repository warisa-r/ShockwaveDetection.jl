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

# This function detects disconuoiities in properties of a flowat a given time step. It computes the gradients of density, velocity, and pressure, 
# and identifies locations where these gradients exceed a threshold (scaled by the maximum gradient of each variable).
function detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x, threshold)

    # Compute first gradients
    density_grad = compute_gradients(density_at_t, x)
    velocity_grad = compute_gradients(velocity_at_t, x)
    pressure_grad = compute_gradients(pressure_at_t, x)

    # Find the maximum gradient for each variable
    max_density_grad = maximum(abs.(density_grad))
    max_velocity_grad = maximum(abs.(velocity_grad))
    max_pressure_grad = maximum(abs.(pressure_grad))

    threshold_density = threshold * max_density_grad
    threshold_velocity = threshold * max_velocity_grad
    threshold_pressure = threshold * max_pressure_grad

    # Find points where the gradient exceeds a certain threshold
    discontinuity_location_density = findall(gradient -> abs(gradient) > threshold_density, density_grad)
    discontinuity_location_velocity = findall(gradient -> abs(gradient) > threshold_velocity, velocity_grad)
    discontinuity_location_pressure = findall(gradient -> abs(gradient) > threshold_pressure, pressure_grad)

    # Combine detections (common shock location across variables)
    shock_locations = intersect(intersect(discontinuity_location_density, discontinuity_location_velocity), discontinuity_location_pressure)
    
    return shock_locations
end

struct GradientShockDetectionAlgo{T} <: Abstract1DShockDetectionAlgo
    threshold::T
end # GradientShockDetectionAlgo

function detect(flow_data::FlowData, alg::GradientShockDetectionAlgo)
    # Unpack all the values from the detector
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field
    bounds = flow_data.bounds
    tsteps = flow_data.tsteps
    threshold = alg.threshold

    len_x = size(density_field, 1)
    x = range(bounds[1][1], stop=bounds[1][2], length=len_x)

    shock_positions_over_time = []
    
    for t_step in 1:length(tsteps)
        density_field_t = density_field[:, t_step]
        velocity_field_t = velocity_field[:, t_step]
        pressure_field_t = pressure_field[:, t_step]
        
        # Use the simple shock detection algorithm to detect the normal shock
        shock_positions = detect_discon_at_timestep(density_field_t, velocity_field_t, pressure_field_t, x, threshold)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end