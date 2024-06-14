# TODO: specify what is used
using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful

"""
    read_output_file(filename)

Reads data from an output file.

# Arguments
- `filename::String`: The path to the output file.

# Returns
- `x0_xmax::Vector{Float64}`: A vector containing the first and last x values.
- `t_values::Vector{Float64}`: A vector containing the time values.
- `u_values::Array{Float64, 3}`: A 3D array containing the u values.
- `dim_u::Vector{Int}`: A vector containing the dimensions of the u values (N_u x N_x).

This function reads data from the specified output file generated from 1d_plots.jl that use Euler2D.jl to solve, which is assumed to have a specific format. 
It reads the dimensions of the u values, the first and last x values, the number of time steps, the time values, and the u values themselves. 
It reshapes the u values into a 3D array and returns all the read data.
"""
function read_output_file(filename)
    open(filename, "r") do f
        dims_u = Vector{Int}(undef, 2)
        read!(f, dims_u)

        # Read the first and last x values
        x0_xmax = Vector{Float64}(undef, 2)
        read!(f, x0_xmax)

        # Read the number of time steps
        num_timesteps = Vector{Int}(undef, 1)
        read!(f, num_timesteps)

        # Read the time values
        t_values = Vector{Float64}(undef, num_timesteps[1])
        read!(f, t_values)

        # Read the u values
        u_values = Vector{Float64}(undef, prod(dims_u)*num_timesteps[1])
        read!(f, u_values)


        # Reshape u_values to a 3D array
        u_values = reshape(u_values, dims_u[1], dims_u[2], num_timesteps[1]) # N_u x N_x x N_t as u a vector of 3 is written in range of x according to each time step


        return x0_xmax, t_values, u_values, dims_u
    end
end

struct FlowData
    u_values::Array{Float64, 3}
    density_field::Array{Float64, 2}
    velocity_field::Array{Float64, 2}
    pressure_field::Array{Float64, 2}
    x0_xmax::Vector{Float64}
    t_values::Array{Float64, 1}
    mach_to_m_s::Bool
end # FlowData

function FlowData(file_path::String, mach_to_m_s=true)
    x0_xmax, t_values, u_values, _ = read_output_file(file_path)
    density_field, velocity_field, pressure_field = convert_to_primitive(u_values, mach_to_m_s)
    return FlowData(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values, mach_to_m_s)
end

# Constructor that takes in the necessary parameters to simulate the 1D flow and returns the FlowData object
# Copy of the function in the script Euler2D.jl/scripts/1d_plots.jl with slight modification
# to be used in the ShockwaveDetection module
# TODO: this doesnt work yet for some reason?
function FlowData(x_min, x_max, ncells_x, x_bcs, T::Float64, u0, 
    gas::CaloricallyPerfectGas = DRY_AIR,
    CFL = 0.75,
    max_tsteps = typemax(Int); verbose = false)

    # Initialize arrays to store results
    u_values = []
    t_values = [0.0]

    # Generate the spatial grid
    xs = range(x_min, x_max; length = ncells_x + 1)
    Δx = step(xs)

    # Initialize the initial condition
    u = stack([u0(x + Δx / 2) for x in xs[1:end-1]])
    u_next = zeros(eltype(u), size(u))

    # Store the initial condition
    push!(u_values, u)

    # Time integration loop
    while ((!(t_values[end] > T || t_values[end] ≈ T)) && length(t_values) <= max_tsteps)
        try
            Δt = maximum_Δt(x_bcs, u, Δx, CFL, 1; gas = gas)
        catch err
            @show length(t_values), t_values[end]
            println("Δt calculation failed.")
            println(typeof(err))
            break
        end

        if t_values[end] + Δt > T
            Δt = T - t_values[end]
        end

        if verbose
            (length(t_values) % 10 == 0) && @show length(t_values), t_values[end], Δt
        end

        # Time step with HLL scheme
        step_euler_hll!(u_next, u, Δt, Δx, x_bcs; gas = gas)

        # Update time and store solution
        push!(t_values, t_values[end] + Δt)
        push!(u_values, u_next)

        # Move to the next time step
        u = u_next
    end

    u_values = cat(u_values..., dims=3)
    x0_xmax = [x_min, x_max]
    density_field, velocity_field, pressure_field = convert_to_primitive(u_values)
    return FlowData(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values)
end