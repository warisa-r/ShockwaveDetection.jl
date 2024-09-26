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
    compare_shock_positions_over_time_1d(original_positions, noisy_positions)

Calculate the mean squared difference between the shock positions from two sets of data over time.

### Parameters:
- `original_positions::Vector{Vector{Int64}}`: A vector of integer vectors representing the original shock positions over time. Each inner vector corresponds to shock positions at a specific time step.
- `noisy_positions::Vector{Vector{Int64}}`: A vector of integer vectors representing the noisy shock positions over time. Each inner vector corresponds to shock positions at the same specific time step.

### Returns:
- `Float64`: The mean squared difference between the shock positions of the two input vectors. If no valid comparisons can be made, the function returns `NaN`.

### Details:
- The function compares each pair of shock position vectors (from `original_positions` and `noisy_positions`). 
- If both vectors for a specific time step are empty, the comparison is skipped.
- If one of the vectors is empty, `Inf` is recorded as a difference.
- For vectors of different lengths, the function calculates the mean squared difference for the overlapping elements.
- After calculating differences, `Inf` values are filtered out before computing the mean of the valid differences.
- The function also tracks how many comparisons resulted in `Inf` and reports the proportion of such comparisons.

### Example:
```julia
original = [[375], [370], [], [367], [366]]
noisy = [[375], [], [368], [367], [365]]

mean_diff = compare_shock_positions_over_time_1d(original, noisy)
println("Mean Squared Difference: ", mean_diff)

"""
function compare_shock_positions_over_time_1d(original_positions, noisy_positions)
    diffs = Float64[]
    inf_count = 0

    println("Length of original_positions: ", length(original_positions))
    println("Length of noisy_positions: ", length(noisy_positions))

    for (orig, noisy) in zip(original_positions, noisy_positions)
        println("Length of current original: ", length(orig))
        println("Length of current noisy: ", length(noisy))

        if isempty(orig) && isempty(noisy)
            println("Both arrays are empty. Skipping comparison.")
            continue
        elseif isempty(orig)
            println("Original array is empty. Recording Inf as a difference.")
            push!(diffs, Inf)
            inf_count += 1
            continue
        elseif isempty(noisy)
            println("Noisy array is empty. Recording Inf as a difference.")
            push!(diffs, Inf)
            inf_count += 1
            continue
        end

        # Calculate the differences for non-empty arrays
        if length(orig) == length(noisy)
            diff = mean((orig .- noisy).^2)
            push!(diffs, diff)
        else
            println("Arrays have different lengths. Original: ", length(orig), " Noisy: ", length(noisy))
            min_len = min(length(orig), length(noisy))
            diff = mean((orig[1:min_len] .- noisy[1:min_len]).^2)
            push!(diffs, diff)
        end
    end

    # Remove any Inf values from diffs before calculating mean
    valid_diffs = filter(!isinf, diffs)

    if isempty(valid_diffs)
        println("No valid differences calculated.")
        return NaN
    end

    mean_diff = mean(valid_diffs)
    inf_ratio = inf_count / length(original_positions)  # Calculate ratio based on original count

    println("Calculated differences: ", diffs)
    println("Mean of valid differences: ", mean_diff)
    println("Number of comparisons with missing data: ", inf_count)
    println("Proportion of comparisons with missing data: ", inf_ratio)

    return mean_diff
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

### Example:
```julia
original_fits = [[Fitting(params1), nothing], [Fitting(params2), Fitting(params3)]]
noisy_fits = [[Fitting(params4)], [Fitting(params5), nothing]]

rmse = compare_shock_fits_over_time(original_fits, noisy_fits)
println("RMSE between fits: ", rmse)
"""
function compare_shock_fits_over_time(original_fit, noisy_fit)
    # Initialize parameter storage
    original_params = Float64[]
    noisy_params = Float64[]
    missing_count = 0

    # Helper function to extract parameters from nested arrays
    function extract_params(fits, params)
        for fit_vector in fits
            for fit in fit_vector
                if fit !== nothing  # Ensure we have a valid Fitting instance
                    push!(params, fit.parameters[1])  # Assuming you want the first parameter
                end
            end
        end
    end

    # Extract parameters for original and noisy fits
    extract_params(original_fit, original_params)
    extract_params(noisy_fit, noisy_params)

    missing_count = abs(length(original_params)-length(noisy_params))
    # Check the extracted parameters
    println("Extracted Original Parameters: ", original_params)
    println("Extracted Noisy Parameters: ", noisy_params)

    println("Length of Original Parameters: ", length(original_params))
    println("Length of Noisy Parameters: ", length(noisy_params))
    println("Number of Missing Parameters: ", missing_count)

    # Ensure both parameter vectors are not empty
    if isempty(original_params) || isempty(noisy_params)
        error("No valid Fitting instances found in the inputs.")
    end

    # Pad the shorter array with zeroes to match lengths
    if length(original_params) > length(noisy_params)
        padding = fill(0.0, length(original_params) - length(noisy_params))
        append!(noisy_params, padding)
    elseif length(noisy_params) > length(original_params)
        padding = fill(0.0, length(noisy_params) - length(original_params))
        append!(original_params, padding)
    end

    # Check the new lengths after padding
    println("Padded Original Parameters: ", original_params)
    println("Padded Noisy Parameters: ", noisy_params)
    println("Length of Padded Original Parameters: ", length(original_params))
    println("Length of Padded Noisy Parameters: ", length(noisy_params))

    # Calculate the differences
    differences = original_params .- noisy_params
    
    # Calculate the RMSE (no NaN values, padded with 0s)
    rmse = sqrt(mean(differences .^ 2))

    return rmse
end