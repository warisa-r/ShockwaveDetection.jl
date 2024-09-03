# TODO: specify what is used
using Euler2D: Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful


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
    noise::Bool
end


function FlowData(file_path::String, mach_to_m_s::Bool=true, noise::Bool=true)
    if endswith(file_path, ".tape")
        sim_data = Euler2D.load_euler_sim(file_path)
        ncells = sim_data.ncells
        nsteps = sim_data.nsteps
        bounds = sim_data.bounds
        tsteps = sim_data.tsteps
        u = sim_data.u
        cell_ids = nothing
        density_field, velocity_field, pressure_field = convert_to_primitive(u, ncells, nsteps, mach_to_m_s, noise)
    elseif endswith(file_path, ".celltape")
        sim_data = Euler2D.load_cell_sim(file_path)
        ncells = sim_data.ncells
        nsteps = sim_data.nsteps
        bounds = sim_data.bounds
        tsteps = sim_data.tsteps
        cell_ids = sim_data.cell_ids
        u = nothing
        density_field, velocity_field, pressure_field = convert_to_primitive(sim_data, nsteps, mach_to_m_s, noise)
    else
        error("Unsupported file type. Please provide a .tape  or .celltape file.")
    end
    return FlowData(ncells, nsteps, bounds, tsteps, u, density_field, velocity_field, pressure_field, mach_to_m_s, cell_ids, noise)
end

# TODO: if the initial condition is given like the scripts, convert EulerSim to FlowData