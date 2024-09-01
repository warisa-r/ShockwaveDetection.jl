using Images
using Statistics: mean

# Function to replace NaNs with the mean of neighboring values for cells near obstacles and obstacles cells
function replace_nan_with_mean!(matrix)
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            if isnan(matrix[i, j])
                # Get the neighborhood (3x3 window) around the NaN value since the kernel size is 3x3
                imin = max(1, i-1)
                imax = min(size(matrix, 1), i+1)
                jmin = max(1, j-1)
                jmax = min(size(matrix, 2), j+1)
                
                # Extract the valid (non-NaN) neighbors
                neighbors = matrix[imin:imax, jmin:jmax]
                valid_neighbors = neighbors[.!isnan.(neighbors)]
                
                # Replace NaN with the mean of valid neighbors
                if !isempty(valid_neighbors)
                    matrix[i, j] = mean(valid_neighbors)
                else # If all neighbors are NaN (in the obstacle zone)
                    matrix[i, j] = 0.0
                end
            end
        end
    end
end


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

function detect_shock_points(flow_data::FlowData, alg::ImageProcessingShockDetectionAlgo, has_obstacle)
    # Unpack all the values from the detector
    nsteps = flow_data.nsteps
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field
    # Preallocate shock_positions_over_time with a size of nsteps, filled with undef
    shock_positions_over_time = Vector{Any}(undef, nsteps)
    
    @threads for t_step in 1:nsteps
        density_field_t = density_field[:, :, t_step]
        velocity_field_x_t = velocity_field[1, :, :, t_step]
        velocity_field_y_t = velocity_field[2, :, :, t_step]
        pressure_field_t = pressure_field[:, :, t_step]

        if has_obstacle
            replace_nan_with_mean!(density_field_t)
            replace_nan_with_mean!(velocity_field_x_t)
            replace_nan_with_mean!(velocity_field_y_t)
            replace_nan_with_mean!(pressure_field_t)
        end

        velocity_field_t = sqrt.(velocity_field_x_t.^2 + velocity_field_y_t.^2) # Calculate magnitude after handling replacement of NaNs
        
        # Use the simple shock detection algorithm to detect the normal shock
        shock_positions_over_time[t_step] = detect_discon_at_timestep(density_field_t, velocity_field_t, pressure_field_t, alg)
    end
    return shock_positions_over_time
end
