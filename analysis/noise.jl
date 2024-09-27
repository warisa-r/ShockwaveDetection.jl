using Euler2D
using Distributions
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
using Statistics  

"""
    struct NoiseData

A structure to store noise configuration for testing and analysis purposes.

# Fields
- `intensity::Float64`: The intensity of the noise.
- `distribution::T`: The probability distribution used to generate noise.
"""
struct NoiseData{T}
    intensity::Float64
    distribution::T
end

"""
    NoiseData(intensity::Float64, distribution::T) -> NoiseData{T}

Constructs a `NoiseData` instance with the specified intensity and distribution.

# Arguments
- `intensity`: A `Float64` representing the intensity of the noise.
- `distribution`: The probability distribution used to generate noise.

# Returns
- A `NoiseData` instance.
"""
function NoiseData(intensity::Float64, distribution::T) where T
    return NoiseData{T}(intensity, distribution)
end

"""
    apply_noise(data::Array{T}, noise_data::NoiseData) -> Array{T}

Applies noise to a given data array using the noise intensity and distribution from `NoiseData`.

# Arguments
- `data`: The original data array to which noise will be applied.
- `noise_data`: A `NoiseData` struct containing the noise configuration.

# Returns
- An array with noise added to the original data.
"""
function apply_noise(data::Array{T}, noise_data::NoiseData) where T
    return data .+ noise_data.intensity * rand(noise_data.distribution, size(data))
end

"""
    apply_noise_to_flow(flow_data::FlowData, noise_data::NoiseData)

Applies noise to the fields of the `FlowData` struct (density, velocity, pressure) using the noise specified in `NoiseData`.

# Arguments
- `flow_data`: The original `FlowData` to which noise will be added.
- `noise_data`: A `NoiseData` struct containing the noise configuration.

# Returns
- A new `FlowData` struct with noise added to the relevant fields.
"""
function apply_noise_to_flow(flow_data::FlowData, noise_data::NoiseData)
    noisy_density_field = apply_noise(flow_data.density_field, noise_data)
    noisy_velocity_field = apply_noise(flow_data.velocity_field, noise_data)
    noisy_pressure_field = apply_noise(flow_data.pressure_field, noise_data)

    # Create a new FlowData struct with noisy fields
    return FlowData(
        flow_data.ncells,
        flow_data.nsteps,
        flow_data.bounds,
        flow_data.tsteps,
        flow_data.u,
        noisy_density_field,
        noisy_velocity_field,
        noisy_pressure_field,
        flow_data.mach_to_m_s,
        flow_data.cell_ids
    )
end

"""
    compare_shock_positions_over_time_1d(original_positions::Vector{Vector{Float64}}, 
                                                       noisy_positions::Vector{Vector{Float64}})

Compare shock positions between original and noisy datasets over time in 1D, using proximity-based matching to better quantify errors caused by false positives, missed detections, or noisy deviations.

# Arguments
- `original_positions::Vector{Vector{Float64}}`: 
    A vector of vectors where each inner vector contains shock positions for a specific time step in the original data. Each time step may have multiple shock positions.
  
- `noisy_positions::Vector{Vector{Float64}}`: 
    A vector of vectors where each inner vector contains shock positions for a specific time step in the noisy data. Similar to `original_positions`, but potentially with noise.

# Returns
- `rmse::Float64`: 
    The Root Mean Squared Error (RMSE) between the matched shock positions from the original and noisy datasets. This metric quantifies the deviation in shock positions due to noise.
  
- `missing_count::Int`: 
    The total count of shock positions that were missing or unmatched between the original and noisy datasets. This includes cases where shocks were present in the original data but absent in the noisy data, or vice versa.

# Methodology
1. For each time step, the function compares shock positions between the `original_positions` and `noisy_positions`.
2. Instead of assuming the shock positions correspond one-to-one in both datasets, each shock in the original data is matched to the **closest** shock in the noisy data. This accounts for variations in the number and order of shock positions due to noise.
3. **Proximity Matching**:
   - The `findmin` function is used to find the closest shock in the noisy data for each shock in the original data. If the distance between the matched positions is within a certain threshold (set to `10` units), the squared difference is computed and stored.
   - If no close match is found, the function applies a **large penalty** (set to `100.0`) to reflect a significant mismatch or false positive.
4. **Handling False Positives**:
   - Any extra shock positions in the noisy data that don't correspond to any shocks in the original data are penalized similarly. These represent **false detections** due to noise.
5. **Handling Missed Detections**:
   - Shock positions in the original data that are missing from the noisy data are also penalized with the large penalty. These represent **missed detections**.
6. After all time steps are processed, the function calculates the RMSE over the matched differences and any penalties.
7. The function also tracks the total count of unmatched (or missing) shocks for additional insight.
"""
# Helper function to find the closest match
function find_closest_match(value, noisy_positions)
    if isempty(noisy_positions)
        return nothing, nothing
    end

    # Find the closest position in noisy_positions
    differences = abs.(noisy_positions .- value)
    min_diff, min_idx = findmin(differences)
    
    # Return the closest noisy position and its index
    return noisy_positions[min_idx], min_idx
end

function compare_shock_positions_over_time_1d(original_positions, noisy_positions, threshold=10.0, large_penalty=100.0)
    diffs = Float64[]
    missing_count = 0

    for (orig, noisy) in zip(original_positions, noisy_positions)
        if isempty(orig) && isempty(noisy)
            # Both are empty, no shocks detected in both original and noisy
            continue
        elseif isempty(noisy) && !isempty(orig)
            # Noisy missed all shocks, apply large penalty for each missing shock
            missing_count += length(orig)
            for _ in orig
                push!(diffs, large_penalty)
            end
            continue
        elseif isempty(orig)
            # No original shock, but noisy detected something (false positives)
            # This may not count as a "missing", but it could be penalized differently
            missing_count += length(noisy)
            for _ in noisy
                push!(diffs, large_penalty)
            end
            continue
        end

        # Match shocks between original and noisy
        for orig_pos in orig
            closest_noisy_pos, closest_idx = find_closest_match(orig_pos, noisy)
            if isnothing(closest_noisy_pos) || abs(orig_pos - closest_noisy_pos) > threshold
                # No match found within threshold, apply penalty
                push!(diffs, large_penalty)
                missing_count += 1
            else
                # Calculate squared difference for valid matches
                diff = (orig_pos - closest_noisy_pos)^2
                push!(diffs, diff)
                # Remove the matched noisy position to avoid double matching
                deleteat!(noisy, closest_idx)
            end
        end
    end

    if isempty(diffs)
        return NaN, missing_count
    end

    rmse = sqrt(mean(diffs))
    return rmse, missing_count
end


"""
    compare_shock_fits_over_time(original_fit, noisy_fit)

Compare the parameters of original and noisy shock fits over time by calculating the root mean square error (RMSE) between them.

### Parameters:
- `original_fit::Vector{Vector{Union{Nothing, Fitting}}}`: A nested vector of fitting instances representing the original shock fits. Each inner vector corresponds to shock fits at a specific time step.
- `noisy_fit::Vector{Vector{Union{Nothing, Fitting}}}`: A nested vector of fitting instances representing the noisy shock fits. Each inner vector corresponds to shock fits at the same specific time step.

### Returns:
- `Float64`: The root mean square error (RMSE) between the extracted parameters of the original and noisy fits.

### Details:
- The function extracts the first parameter from each fitting instance (if valid) and compares the parameters from both original and noisy fits.
- If there are any missing parameters (instances of `nothing`), they are counted and reported.
- If the lengths of the parameter vectors differ, the shorter vector is padded with zeroes to match the length of the longer vector.
- The function computes the RMSE based on the differences between the parameter vectors.
"""
function compare_shock_fits_over_time(original_fit, noisy_fit)
    # Initialize parameter storage and missing count
    original_params = Float64[]
    noisy_params = Float64[]
    missing_count = 0  # Count how many parameters are missing

    # Helper function to extract parameters from nested arrays
    function extract_params(fits, params)
        for fit_vector in fits
            for fit in fit_vector
                if fit !== nothing
                    # Append all elements of fit.parameters to params
                    append!(params, fit.parameters)  
                end
            end
        end
    end

    # Extract parameters for original and noisy fits
    extract_params(original_fit, original_params)
    extract_params(noisy_fit, noisy_params)

    # Check the extracted parameters
    println("Extracted Original Parameters: ", original_params)
    println("Extracted Noisy Parameters: ", noisy_params)

    # Find the maximum length of parameter lists (for padding)
    max_len = max(length(original_params), length(noisy_params))

    # Pad the shorter array with NaNs (or another value like Inf) to indicate missing parameters
    if length(original_params) < max_len
        missing_count += max_len - length(original_params)  # Increment missing count
        padding = fill(NaN, max_len - length(original_params))
        append!(original_params, padding)
    end

    if length(noisy_params) < max_len
        missing_count += max_len - length(noisy_params)  # Increment missing count
        padding = fill(NaN, max_len - length(noisy_params))
        append!(noisy_params, padding)
    end

    # Check the new lengths after padding
    println("Padded Original Parameters: ", original_params)
    println("Padded Noisy Parameters: ", noisy_params)

    # Calculate the differences while ignoring NaNs
    differences = original_params .- noisy_params
    filtered_differences = filter(!isnan, differences)  # Remove NaNs for RMSE calculation

    if isempty(filtered_differences)
        println("No valid parameter differences to calculate RMSE.")
        return NaN, missing_count  # If all values were missing
    end

    # Calculate the RMSE (ignoring NaN values)
    rmse = sqrt(mean(filtered_differences .^ 2))

    return rmse, missing_count
end
