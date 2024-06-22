using Images

struct ImageProcessingShockDetectionAlgo <: Abstract2DShockDetectionAlgo
    threshold::Float64
    kernelname::Symbol
end # ImageProcessingShockDetectionAlgo

function get_kernel_by_name(kernelname::Symbol)
    if kernelname == :sobel
        return KernelFactors.sobel
    elseif kernelname == :prewitt
        return KernelFactors.prewitt
    elseif kernelname == :ando3
        return KernelFactors.ando3
    elseif kernelname == :scharr
        return KernelFactors.scharr
    elseif kernelname == :bickley
        return KernelFactors.bickley
    else
        error("Unsupported kernel name: $kernelname")
    end
end

# Function to detect sharp gradients
function detect_discon_at_timestep(density_at_t, velocity_at_t, pressure_at_t, alg::ImageProcessingShockDetectionAlgo)
    threshold = alg.threshold
    kernelfunc = get_kernel_by_name(alg.kernelname)

    # Compute first gaussian derivatives along x and y directions for all property fields
    density_grad_dx, density_grad_dy = imgradients(density_at_t, kernelfunc, "replicate")
    velocity_grad_dx, velocity_grad_dy = imgradients(velocity_at_t, kernelfunc, "replicate")
    pressure_grad_dx, pressure_grad_dy = imgradients(pressure_at_t, kernelfunc, "replicate")

    # Calculate gradient magnitudes for each property
    density_gradient_magnitudes = sqrt.(density_grad_dx.^2 + density_grad_dy.^2)
    velocity_gradient_magnitudes = sqrt.(velocity_grad_dx.^2 + velocity_grad_dy.^2)
    pressure_gradient_magnitudes = sqrt.(pressure_grad_dx.^2 + pressure_grad_dy.^2)

    # Find maximum gradient magnitudes
    max_density_gradient_magnitude = maximum(density_gradient_magnitudes)
    max_velocity_gradient_magnitude = maximum(velocity_gradient_magnitudes)
    max_pressure_gradient_magnitude = maximum(pressure_gradient_magnitudes)

    # Scale the maximum values to use as thresholds
    density_threshold = max_density_gradient_magnitude * threshold
    velocity_threshold = max_velocity_gradient_magnitude * threshold
    pressure_threshold = max_pressure_gradient_magnitude * threshold

    # Apply threshold to detect sharp gradients for each property
    density_sharp_gradients = findall(gradient -> abs(gradient) > density_threshold, density_gradient_magnitudes)
    velocity_sharp_gradients = findall(gradient -> abs(gradient) > velocity_threshold, velocity_gradient_magnitudes)
    pressure_sharp_gradients = findall(gradient -> abs(gradient) > pressure_threshold, pressure_gradient_magnitudes)

    # Assuming the variables are arrays or sets of indices
    high_gradient_intersection = intersect(density_sharp_gradients, velocity_sharp_gradients, pressure_sharp_gradients)
    
    # return n_shock-element Vector{CartesianIndex{2}}
    return high_gradient_intersection
end

function detect(flow_data::FlowData, alg::ImageProcessingShockDetectionAlgo)
    # Unpack all the values from the detector
    nsteps = flow_data.nsteps
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field
    shock_positions_over_time = []
    
    for t_step in 1:nsteps
        density_field_t = density_field[:, :, t_step]
        velocity_field_x_t = velocity_field[1, :, :, t_step]
        velocity_field_y_t = velocity_field[2, :, :, t_step]
        velocity_field_t = sqrt.(velocity_field_x_t.^2 + velocity_field_y_t.^2)
        pressure_field_t = pressure_field[:, :, t_step]
        
        # Use the simple shock detection algorithm to detect the normal shock
        shock_positions = detect_discon_at_timestep(density_field_t, velocity_field_t, pressure_field_t, alg)
        push!(shock_positions_over_time, shock_positions)
    end
    return shock_positions_over_time
end
