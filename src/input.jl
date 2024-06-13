"""
    read_output_file(filename)

Reads data from an output file.

# Arguments
- `filename::String`: The path to the output file.

# Returns
- `x0_xmax::Vector{Float64}`: A vector containing the first and last x values.
- `t_values::Vector{Float64}`: A vector containing the time values.
- `u_values::Array{Float64, 3}`: A 3D array containing the u values.
- `dims_u::Vector{Int}`: A vector containing the dimensions of the u values (N_u x N_x).

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
end # FlowData

function FlowData(file_path::String)
    x0_xmax, t_values, u_values, _ = read_output_file(file_path)
    density_field, velocity_field, pressure_field = convert_to_primitive(u_values)
    return FlowData(u_values, density_field, velocity_field, pressure_field, x0_xmax, t_values)
end

#TODO: Implement this
# Call the solver with given boundary conditions
function FlowData()
    println("Not implemented yet")
end