# Use Interpolations.jl to perform np.gradient-like operations
using LinearAlgebra: I
using Unitful: ustrip
using Interpolations: interpolate, Gridded, Linear, gradient
using ShockwaveProperties: state_behind, ConservedProps, DRY_AIR, mach_number, shock_temperature_ratio, shock_density_ratio, shock_pressure_ratio

abstract type AbstractShockDetectionAlgo end
# Abstract base type for 1D shock detectors
abstract type Abstract1DShockDetectionAlgo <: AbstractShockDetectionAlgo end

# Function to compute gradients using Interpolations.jl
# without explicit function definition
function compute_gradients(arr::AbstractVector, x)
    # Convert StepRangeLen to LinRange because that's what Interpolations.jl expects
    x_linrange = LinRange(first(x), last(x), length(x))
    itp = interpolate((x_linrange,), arr, Gridded(Linear()))
    grad = only.(gradient.(Ref(itp), x_linrange))
    return grad
end

#TODO: more efficient?
function shock_rankine_hugoniot(shock_locations, u_values_t)
    # Define normal and tangential vectors for 1D case because state_behind needs it
    n = [1]
    t = [0]

    normal_shock_positions = []

    # flux for the euler equations from ShockwaveProperties.jl's test since it's not exported
    function F(u::ConservedProps)
        v = u.ρv / u.ρ # velocity
        P = pressure(u; gas = DRY_AIR)
        return vcat(u.ρv', (u.ρv .* v' + I * P), (v .* (u.ρE + P))') # stack row vectors
    end

    for shock_location in shock_locations
        u_L = ConservedProps(u_values_t[:, shock_location])
        # Check if Mach number is greater than 1 at the shock location -> for normal shock in supersonic flow
        # else, the shock is not a normal shock and this is not considered
        if mach_number(u_L; gas = DRY_AIR)[1] >= 1
            # TODO: Look further away to get left and right state
            u_R = state_behind(u_L, n, t; gas = DRY_AIR)
            # Check Rankine-Hugonoit conditions
            if all(isapprox.(ustrip.(F(u_L) .* n), ustrip.(F(u_R) .* n); atol=1e-6))
                println("Rankine-Hugonoit conditions are satisfied at shock location $shock_location")
                push!(normal_shock_positions, shock_location)
            else
                # Debug message
                println("Rankine-Hugonoit conditions are not satisfied at shock location $shock_location, mach number = $(mach_number(u_L; gas = DRY_AIR))")
            end
        else
            println("Mach number is less than 1")
        end

        
        
    end

    return normal_shock_positions

end

# This function detects shock locations in a fluid at a given time step. It computes the gradients of density, velocity, and pressure, 
# and identifies locations where these gradients exceed a threshold (scaled by the maximum gradient of each variable). 
# It then checks if the Rankine-Hugoniot conditions are satisfied at these locations.
function detect_normal_shocks_at_timestep(u_values_t, density_at_t, velocity_at_t, pressure_at_t, x, threshold)

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
    shock_location_density = findall(gradient -> abs(gradient) > threshold_density, density_grad)
    shock_location_velocity = findall(gradient -> abs(gradient) > threshold_velocity, velocity_grad)
    shock_location_pressure = findall(gradient -> abs(gradient) > threshold_pressure, pressure_grad)

    # Combine detections (common shock location across variables)
    shock_locations = intersect(intersect(shock_location_density, shock_location_velocity), shock_location_pressure)

    #  Check if Rankine-Hugoniot conditions are satisfied
    shock_locations = shock_rankine_hugoniot(shock_locations, u_values_t)
    
    return shock_locations
end

# TODO: test if the version with the constructor work, if it does then this is probably redundant
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

struct Simple1DShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    u_values::Array{Float64, 3}
    density_field::Array{Float64, 2}
    velocity_field::Array{Float64, 2}
    pressure_field::Array{Float64, 2}
    x0_xmax::Vector{Float64}
    t_values::Array{Float64, 1}
end # Simple1DShockDetector

# Constructor to read values of u in a grid from the output file generated by Euler2d.jl/scripts/1D_plots.jl
function Simple1DShockDetectionAlgo(file_path::String)
    x0_xmax, t_values, u_values, _ = read_output_file(file_path)
    density_field, velocity_field, pressure_field = convert_to_primitive(u_values)
    return Simple1DShockDetectionAlgo(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values)
end

function detect(alg::Simple1DShockDetectionAlgo, threshold = 0.5)
    # Unpack all the values from the detector
    density_field = alg.density_field
    velocity_field = alg.velocity_field
    pressure_field = alg.pressure_field
    x0_xmax = alg.x0_xmax
    t_values = alg.t_values
    u_values = alg.u_values

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