using Euler2D
using Distributions
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

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

Compares shock positions over time between original and noisy data. The function checks if both arrays are empty and skips comparisons accordingly. If one array is empty, it also skips that comparison.
When the lengths differ, it computes the mean difference only for the overlapping elements.

# Arguments
- `original_positions`: Shock positions from the original `FlowData`.
- `noisy_positions`: Shock positions from the noisy `FlowData`.

# Returns
- A measure of the difference between the original and noisy positions.
"""
function compare_shock_positions_over_time_1d(original_positions, noisy_positions)
    diffs = Float64[]

    println("Length of original_positions: ", length(original_positions))
    println("Length of noisy_positions: ", length(noisy_positions))

    for (orig, noisy) in zip(original_positions, noisy_positions)
        println("Length of current original: ", length(orig))
        println("Length of current noisy: ", length(noisy))

        if isempty(orig) && isempty(noisy)
            println("Both arrays are empty. Skipping comparison.")
            continue
        elseif isempty(orig) || isempty(noisy)
            println("One of the arrays is empty. Skipping comparison. Original: ", orig, " Noisy: ", noisy)
            continue
        end

        # Calculate the differences
        if length(orig) == length(noisy)
            diff = mean((orig .- noisy).^2)
            push!(diffs, diff)
        else
            println("Arrays have different lengths. Original: ", length(orig), " Noisy: ", length(noisy))
            # Optionally handle this case, e.g., calculate differences for the overlapping indices
            min_len = min(length(orig), length(noisy))
            diff = mean((orig[1:min_len] .- noisy[1:min_len]).^2)
            push!(diffs, diff)
        end
    end

    if isempty(diffs)
        println("No valid differences calculated.")
        return NaN
    end

    println("Calculated differences: ", diffs)

    return mean(diffs)
end

"""
    compare_shock_positions_over_time_2d(original_positions, noisy_positions)

Compares shock positions over time between original and noisy data. The function checks if both arrays are empty and skips comparisons accordingly. If one array is empty, it also skips that comparison.
When the lengths differ, it computes the mean difference only for the overlapping elements. Additionally, the Tuple(ci) converts a CartesianIndex object like CartesianIndex into a tuple, allowing for easy element-wise operations.

# Arguments
- `original_positions`: Shock positions from the original `FlowData`.
- `noisy_positions`: Shock positions from the noisy `FlowData`.

# Returns
- A measure of the difference between the original and noisy positions.
"""
function compare_shock_positions_over_time_2d(original_positions, noisy_positions)
    diffs = Float64[]

    println("Length of original_positions: ", length(original_positions))
    println("Length of noisy_positions: ", length(noisy_positions))

    for (orig, noisy) in zip(original_positions, noisy_positions)
        println("Length of current original: ", length(orig))
        println("Length of current noisy: ", length(noisy))

        if isempty(orig) && isempty(noisy)
            println("Both arrays are empty. Skipping comparison.")
            continue
        elseif isempty(orig) || isempty(noisy)
            println("One of the arrays is empty. Skipping comparison. Original: ", orig, " Noisy: ", noisy)
            continue
        end

        # Convert CartesianIndex objects to tuples
        orig_positions = [Tuple(ci) for ci in orig]
        noisy_positions = [Tuple(ci) for ci in noisy]

        # Calculate the differences
        if length(orig_positions) == length(noisy_positions)
            diff = mean(sum((orig_pos .- noisy_pos).^2) for (orig_pos, noisy_pos) in zip(orig_positions, noisy_positions))
            push!(diffs, diff)
        else
            println("Arrays have different lengths. Original: ", length(orig_positions), " Noisy: ", length(noisy_positions))
            min_len = min(length(orig_positions), length(noisy_positions))
            diff = mean(sum((orig_positions[i] .- noisy_positions[i]).^2) for i in 1:min_len)
            push!(diffs, diff)
        end
    end

    if isempty(diffs)
        println("No valid differences calculated.")
        return NaN
    end

    println("Calculated differences: ", diffs)

    return mean(diffs)
end



"""
    compare_shock_clusters_over_time(original_clusters, noisy_clusters)

Compares shock clusters over time between original and noisy data.

# Arguments
- `original_clusters`: Shock clusters from the original `FlowData`.
- `noisy_clusters`: Shock clusters from the noisy `FlowData`.

# Returns
- A measure of the difference between the original and noisy clusters.
"""
function compare_shock_clusters_over_time(original_clusters, noisy_clusters)
    diff = original_clusters .- noisy_clusters
    return sqrt(mean(diff.^2))
end

"""
    compare_shock_fits_over_time(original_fits, noisy_fits)

Compares shock fits over time between original and noisy data.

# Arguments
- `original_fits`: Shock fits from the original `FlowData`.
- `noisy_fits`: Shock fits from the noisy `FlowData`.

# Returns
- A measure of the difference between the original and noisy fits.
"""
function compare_shock_fits_over_time(original_fits, noisy_fits)
    diff = original_fits .- noisy_fits
    return sqrt(mean(diff.^2))
end

