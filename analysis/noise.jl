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

# Helper function to find the closest match for a shock in the noisy data
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

# Main function to compare 1D shock positions over time
function compare_shock_positions_over_time_1d(original_positions, noisy_positions; threshold=10.0)
    diffs = Float64[]  # To store the differences between shocks
    missing_count = 0  # Track the number of missing shocks

    for (orig, noisy) in zip(original_positions, noisy_positions)
        if isempty(orig) && isempty(noisy)
            # No shocks in both original and noisy, no comparison needed
            continue
        elseif isempty(noisy) && !isempty(orig)
            # No shocks detected in noisy but there are shocks in original
            missing_count += length(orig)  # Count all as missing
            continue
        elseif isempty(orig) && !isempty(noisy)
            # Original has no shocks, but noisy has extra shocks (false positives)
            missing_count += length(noisy)  # Count all noisy as false detections
            continue
        end

        # Match shocks between original and noisy
        for orig_pos in orig
            closest_noisy_pos, closest_idx = find_closest_match(orig_pos, noisy)

            if isnothing(closest_noisy_pos) || abs(orig_pos - closest_noisy_pos) > threshold
                # No match found within the threshold, consider it missing
                missing_count += 1
            else
                # Calculate squared difference for valid matches
                diff = (orig_pos - closest_noisy_pos)^2
                push!(diffs, diff)
                # Remove the matched noisy position to avoid double matching
                deleteat!(noisy, closest_idx)
            end
        end

        # Any remaining noisy positions are unmatched, count as false positives
        if !isempty(noisy)
            missing_count += length(noisy)
        end
    end

    # If there are no valid differences to compute RMSE, return NaN and the count of missing positions
    if isempty(diffs)
        return NaN, missing_count
    end

    # Calculate RMSE only on the valid matched differences
    rmse = sqrt(mean(diffs))
    return rmse, missing_count
end

"""
"""
# Helper function to calculate residuals between fit and original data points from the shock point detection
function calculate_residuals(fit::ShockwaveDetection.Fitting, data_points)
    residuals = Float64[]

    for i in eachindex(data_points) 
        fitted_value = fit.model(data_points[i][1], fit.parameters)  # Call the model function with parameters
        actual_value = data_points[i][2]                             # Actual y-value from the data
        println("Fitted value: ", fitted_value, " Actual value: ", actual_value)
        push!(residuals, abs(fitted_value - actual_value))
    end

    return residuals
end

# Main function to compare fits over time by analyzing residuals rather than raw parameters
function compare_shock_fits_over_time(original_fit, noisy_fit, original_shock_positions_over_time, noisy_shock_positions_over_time, flow_data; threshold=10.0)
    bounds = flow_data.bounds
    ncells = flow_data.ncells
    nsteps = flow_data.nsteps
    rmse_results = Vector{Float64}(undef, nsteps)  
    data_points = Vector{Any}(undef, nsteps)       

    # Define the x-range
    x = range(bounds[1][1], bounds[1][2], length=ncells[1])

    # Define the y-range
    y = range(bounds[2][1], bounds[2][2], length=ncells[2])

    for t in eachindex(nsteps)  # Use eachindex for time steps
        if isempty(original_shock_positions_over_time[t]) 
            rmse_results[t] = NaN  # No shock positions, set RMSE as NaN
            data_points[t] = []     # No data points to analyze
            continue
        end

        # Transform shock positions to Cartesian coordinates
        original_data_points[t] = ShockwaveDetection.cartesian_index_to_xy(original_shock_positions_over_time[t], x, y)
        noisy_data_points[t] = ShockwaveDetection.cartesian_index_to_xy(noisy_shock_positions_over_time[t], x, y)

        # Calculate residuals for original and noisy fits
        original_residuals = calculate_residuals(original_fit[t], original_data_points[t])  
        noisy_residuals = calculate_residuals(noisy_fit[t], noisy_data_points[t])       

        if isempty(original_residuals) || isempty(noisy_residuals)
            error("No valid residuals found for the input fits.")
        end

        # Compare the residuals to see how well the fits match the actual data points
        residual_diffs = abs.(original_residuals .- noisy_residuals)

        # Calculate the root mean square error between residuals
        rmse_results[t] = sqrt(mean(residual_diffs.^2))
    end

    return rmse_results
end