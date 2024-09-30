using Euler2D
using Distributions
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection
using Statistics
using .NoiseAnalysis

"""
    struct NoiseData

A structure to store noise configuration for testing and analysis purposes.

# Fields
- `intensity::T`: The intensity of the noise.
- `distribution::D`: The probability distribution used to generate noise.
"""
struct NoiseData{T, D}
    intensity::T
    distribution::D
end

"""
    NoiseData(intensity::T, distribution::D) -> NoiseData{T}

Constructs a `NoiseData` instance with the specified intensity and distribution.

# Arguments
- `intensity`: Representing the intensity of the noise.
- `distribution`: The probability distribution used to generate noise.

# Returns
- A `NoiseData` instance.
"""
function NoiseData(intensity::T, distribution::D) where {T, D}
    return NoiseData{T, D}(intensity, distribution)
end

"""
    apply_noise(data::Array{T}, noise_data::NoiseData) -> Array{T, D}

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
    find_closest_match(value, noisy_positions)

Helper function to find the closest match to a given `value` in `noisy_positions`.

# Arguments
- `value`: A numeric value representing the original position.
- `noisy_positions`: A vector of numeric values representing the noisy positions.

# Returns
- The closest matching value from `noisy_positions` and its index.
- Returns `nothing, nothing` if `noisy_positions` is empty.
"""
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

"""
    compare_shock_positions_over_time_1d(original_positions, noisy_positions; threshold::T=10.0)

Main function to compare 1D shock positions between original and noisy data over time.
It attempts to match shocks in `original_positions` to those in `noisy_positions` and computes the RMSE (Root Mean Square Error) between matched shocks.
Unmatched shocks are counted as either false positives (extra detections) or false negatives (missed detections).

# Arguments
- `original_positions`: A vector of vectors containing shock positions from the original data for each time step.
- `noisy_positions`: A vector of vectors containing shock positions from noisy data for each time step.
- `threshold`: A numeric value (default 10.0) specifying the maximum allowable difference between matched shocks.

# Returns
- `rmse`: The Root Mean Square Error between matched shock positions, or `NaN` if no valid matches exist.
- `false_negatives_count`: The number of shocks in the original data that were missed in the noisy data.
- `false_positives_count`: The number of extra shocks detected in the noisy data that do not correspond to any shocks in the original data.
"""
function compare_shock_positions_over_time_1d(original_positions, noisy_positions; threshold::T = 10.0) where T
    diffs = Any[]  # To store the differences between matched shocks
    false_negatives_count = 0  # Count of missed shocks (false negatives)
    false_positives_count = 0  # Count of extra shocks (false positives)

    for (orig, noisy) in zip(original_positions, noisy_positions)
        if isempty(orig) && isempty(noisy)
            # No shocks in both original and noisy, no comparison needed
            continue
        elseif isempty(noisy) && !isempty(orig)
            # No shocks detected in noisy but there are shocks in original (false negatives)
            false_negatives_count += length(orig)  # Count all as missed
            continue
        elseif isempty(orig) && !isempty(noisy)
            # Original has no shocks, but noisy has extra shocks (false positives)
            false_positives_count += length(noisy)  # Count all noisy as false detections
            continue
        end

        # Match shocks between original and noisy
        for orig_pos in orig
            closest_noisy_pos, closest_idx = find_closest_match(orig_pos, noisy)

            if isnothing(closest_noisy_pos) || abs(orig_pos - closest_noisy_pos) > threshold
                # No match found within the threshold, count as false negative
                false_negatives_count += 1
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
            false_positives_count += length(noisy)
        end
    end

    # If there are no valid differences to compute RMSE, return NaN and the counts
    if isempty(diffs)
        return NaN, false_negatives_count, false_positives_count
    end

    # Calculate RMSE only on the valid matched differences
    rmse = sqrt(mean(diffs))
    return rmse, false_negatives_count, false_positives_count
end

"""
    calculate_residuals(fit::ShockwaveDetection.Fitting, data_points)

Calculates the residuals between the actual data points and the fitted values generated by the model.

# Arguments
- `fit`: A `ShockwaveDetection.Fitting` object that contains the model and its parameters used to generate fitted values.
- `data_points`: A vector of tuples or arrays, where each element is a pair `[x_value, actual_value]`. The `x_value` represents the independent variable, and `actual_value` is the observed data point.

# Returns
- `residuals`: A vector of residuals, where each residual is the absolute difference between the fitted value from the model and the actual observed value.

# Method
1. For each data point, extract the `x_value` (independent variable) and `actual_value` (dependent variable).
2. Pass the `x_value` to the model using the parameters from `fit`.
3. Calculate the residual as the absolute difference between the model's fitted value and the `actual_value`.
4. Return the vector of residuals.
"""
function calculate_residuals(fit::ShockwaveDetection.Fitting{T}, data_points::Vector{Vector{T}}) where T
    residuals = Any[]

    for point in data_points
        x_value = point[1]  
        actual_value = point[2]  

        # Call the model function with parameters
        fitted_value_vector = fit.model([x_value], fit.parameters)  
        fitted_value = fitted_value_vector[1]  

        #println("Fitted value: ", fitted_value, " Actual value: ", actual_value)

        # Calculate residual
        push!(residuals, abs(fitted_value - actual_value))  # Calculate residual
    end

    return residuals
end

"""
    compare_shock_fits_over_time(original_fit, noisy_fit, original_shock_positions_over_time, noisy_shock_positions_over_time, flow_data; threshold=10.0)

Compares fitted curves (shock fits) to the respective clusters of shock points over time.

# Arguments
- `original_fit`: A vector of fitted curves from the original data over time.
- `noisy_fit`: A vector of fitted curves from the noisy data over time.
- `original_shock_positions_over_time`: A vector of shock positions from the original data for each time step.
- `noisy_shock_positions_over_time`: A vector of shock positions from the noisy data for each time step.
- `flow_data`: An object containing properties of flow data such as bounds, number of cells, and number of time steps.
- `threshold`: A numeric value (default 10.0) specifying the maximum allowable difference between matched shocks.

# Returns
- `rmse_results`: A vector containing the root mean square error (RMSE) or norm differences of the fit comparisons over time.
- `missing_counts`: A dictionary containing counts of mismatches in different categories:
  - `fit_mismatch_count`: Number of mismatches between the number of fits and clusters.
  - `residual_mismatch_count`: Number of mismatches between the number of residuals calculated.
  - `norm_mismatch_count`: Number of mismatches resolved by comparing the second norm of residuals.

# Method
1. For each time step, the function iterates over the fitted curves and their respective clusters of shock points.
2. It calculates the residuals for both the original and noisy fits.
3. If the residuals differ in length, it calculates the difference in their L2 norms instead of performing RMSE.
4. It assesses the matches based on the specified `threshold` and counts mismatches in various categories.
5. The function compiles and returns the RMSE or norm differences for each time step along with counts of different types of mismatches for further analysis.
"""

function compare_shock_fits_over_time(
    original_result::ShockwaveDetection.ShockDetectionResult2D, 
    noisy_result::ShockwaveDetection.ShockDetectionResult2D, 
    flow_data; threshold=10.0) 

    bounds = flow_data.bounds
    nsteps = flow_data.nsteps

    # Preallocate RMSE results for each time step
    rmse_results = Vector{Any}(undef, nsteps)
    data_points = Vector{Any}(undef, nsteps)

    # Counters for different types of mismatches
    cluster_fit_mismatch_count = 0          # Number of times clusters don't match fits
    empty_residuals_count = 0               # Number of times residuals were empty
    shock_data_missing_count = 0            # Number of times shock points were missing
    norm_mismatch_count = 0                 # Number mismatches resolved by comparing the second norm of residuals.

    # Use the shock clusters from the ShockDetectionResult2D struct
    original_clusters = original_result.shock_clusters_over_time
    noisy_clusters = noisy_result.shock_clusters_over_time

    # Check for mismatch in the number of total clusters and fits before time step iteration
    if length(original_clusters) != length(original_result.shock_fits_over_time)
        error("Mismatch between number of clusters and fits in original data.")
    end
    if length(noisy_clusters) != length(noisy_result.shock_fits_over_time)
        error("Mismatch between number of clusters and fits in noisy data.")
    end

    for t in 1:nsteps
        # Handle cases where shock data is missing in one or both datasets
        if isempty(original_result.shock_positions_over_time[t]) || isempty(noisy_result.shock_positions_over_time[t])
            shock_data_missing_count += 1
            rmse_results[t] = NaN
            data_points[t] = []
            continue
        end

        num_original_clusters = length(original_clusters[t])
        num_noisy_clusters = length(noisy_clusters[t])
        num_fits = length(original_result.shock_fits_over_time[t])  # Number of fits at time step t

        # Check if the number of clusters matches the number of fits for both datasets
        if num_fits != num_original_clusters || num_fits != num_noisy_clusters
            cluster_fit_mismatch_count += 1
            rmse_results[t] = NaN
            data_points[t] = []
            continue
        end

        for i in 1:num_fits
            # Extract shock cluster data for the current fit
            original_data_points = original_clusters[t][i]
            noisy_data_points = noisy_clusters[t][i]

            # Calculate residuals for original and noisy fits
            original_residuals = calculate_residuals(original_result.shock_fits_over_time[t][i], original_data_points)
            noisy_residuals = calculate_residuals(noisy_result.shock_fits_over_time[t][i], original_data_points)

            # If residuals are empty, increment the empty residuals counter
            if isempty(original_residuals) || isempty(noisy_residuals)
                empty_residuals_count += 1
                continue
            end

            # Ensure residuals have the same length
            if length(original_residuals) != length(noisy_residuals)
                # Calculate the second norms of both residuals
                original_norm = norm(original_residuals, 2)
                noisy_norm = norm(noisy_residuals, 2)

                # Calculate the difference between the norms
                norm_difference = abs(original_norm - noisy_norm)

                # Store the difference instead of RMSE
                rmse_results[t] = norm_difference
                data_points[t] = []

                norm_mismatch_count += 1  
                continue
            end

            # Compare the residuals between original and noisy fits
            residual_diffs = abs.(original_residuals .- noisy_residuals)

            # Calculate root mean square error (RMSE)
            rmse_results[t] = sqrt(mean(residual_diffs.^2))
        end
    end

    # Return the RMSE results along with detailed mismatch counts
    return rmse_results, cluster_fit_mismatch_count, empty_residuals_count, shock_data_missing_count, norm_mismatch_count
end