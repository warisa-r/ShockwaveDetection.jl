# TODO: specify what is used
using Euler2D: Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful

struct FlowData{N,NAXES,T}
    ncells::NTuple{N,Int}
    nsteps::Int
    bounds::NTuple{N,Tuple{T,T}}
    tsteps::Vector{T}
    u::Array{T,NAXES}
    density_field::Matrix{T}
    velocity_field::Matrix{T}
    pressure_field::Matrix{T}
    mach_to_m_s::Bool
end

function FlowData(file_path::String, mach_to_m_s=true)
    if endswith(file_path, ".tape")
        sim_data = Euler2D.load_euler_sim(file_path)
        ncells = sim_data.ncells
        nsteps = sim_data.nsteps
        bounds = sim_data.bounds
        tsteps = sim_data.tsteps
        u = sim_data.u
        density_field, velocity_field, pressure_field = convert_to_primitive(u, mach_to_m_s)
    else
        error("Unsupported file type. Please provide a .out or .tape file.")
    end
    return FlowData(ncells, nsteps, bounds, tsteps, u, density_field, velocity_field, pressure_field, mach_to_m_s)
end

# TODO: if the initial condition is given like the scripts, convert EulerSim to FlowData