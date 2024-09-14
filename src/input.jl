# TODO: specify what is used
using Euler2D: Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful

"""
    struct FlowData{N, T}

A structure to hold flow data for simulations.

# Fields
- `ncells::NTuple{N, Int}`: Number of cells in each dimension.
- `nsteps::Int`: Number of time steps in the simulation.
- `bounds::NTuple{N, Tuple{T, T}}`: Bounds of the simulation domain in each dimension.
- `tsteps::Vector{T}`: Time steps of the simulation.
- `u::Union{Array{T}, Nothing}`: Array of primitive variables or `Nothing` for `.celltape` data.
- `density_field::Array{T}`: Density field of the simulation.
- `velocity_field::Array{T}`: Velocity field of the simulation.
- `pressure_field::Array{T}`: Pressure field of the simulation.
- `mach_to_m_s::Bool`: Flag indicating if Mach number should be converted to meters per second.
- `cell_ids::Union{Matrix{Int64}, Nothing}`: Cell IDs for `.celltape` data. `Nothing` for `.tape` data.
"""
struct FlowData{N, T}
    ncells::NTuple{N, Int}
    nsteps::Int
    bounds::NTuple{N, Tuple{T, T}}
    tsteps::Vector{T}
    u::Union{Array{T}, Nothing} # Nothing for .celltape data since all other important data are stored in other attributes already
    density_field::Array{T}
    velocity_field::Array{T}
    pressure_field::Array{T}
    mach_to_m_s::Bool
    cell_ids::Union{Matrix{Int64}, Nothing}
end

"""
    FlowData(file_path::String, mach_to_m_s=true)

Load flow data from a file and create a `FlowData` object.

# Arguments
- `file_path::String`: Path to the file containing the simulation data.
- `mach_to_m_s::Bool=true`: Flag indicating if Mach number should be converted to meters per second.

# Returns
- `FlowData`: A `FlowData` object containing the loaded simulation data.

# Description
This function loads flow data from a specified file path. It supports both `.tape` and `.celltape` file formats. For `.tape` files, it loads Euler simulation data and converts it to primitive variables. For `.celltape` files, it loads cell simulation data.
"""
function FlowData(file_path::String, mach_to_m_s=true)
    if endswith(file_path, ".tape")
        sim_data = Euler2D.load_euler_sim(file_path)
        ncells = sim_data.ncells
        nsteps = sim_data.nsteps
        bounds = sim_data.bounds
        tsteps = sim_data.tsteps
        u = sim_data.u
        cell_ids = nothing
        density_field, velocity_field, pressure_field = convert_to_primitive(u, ncells, nsteps, mach_to_m_s)
    elseif endswith(file_path, ".celltape")
        # for celltape files, velocity field is currently already m_s. No feature to convert it to mach number.
        sim_data = Euler2D.load_cell_sim(file_path)
        ncells = sim_data.ncells
        nsteps = sim_data.nsteps
        bounds = sim_data.bounds
        tsteps = sim_data.tsteps
        cell_ids = sim_data.cell_ids
        u = nothing
        density_field, velocity_field, pressure_field = convert_to_primitive(sim_data, nsteps)
    else
        error("Unsupported file type. Please provide a .tape  or .celltape file.")
    end
    return FlowData(ncells, nsteps, bounds, tsteps, u, density_field, velocity_field, pressure_field, mach_to_m_s, cell_ids)
end

# TODO: if the initial condition is given like the scripts, convert EulerSim to FlowData