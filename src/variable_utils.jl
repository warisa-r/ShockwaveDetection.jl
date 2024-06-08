using ShockwaveProperties: primitive_state_vector, pressure, speed_of_sound, DRY_AIR
using Unitful: ustrip

"""
    convert_to_primitive(u_values)

Converts conserved state variables to primitive variables.

# Arguments
- `u_values::Array{Float64, 3}`: A 3D array containing the conserved state variables (density, momentum, total energy) at each point x and time t.

# Returns
- `density_field::Array{Float64, 2}`: A 2D array containing the density field.
- `velocity_field::Array{Float64, 2}`: A 2D array containing the velocity field.
- `pressure_field::Array{Float64, 2}`: A 2D array containing the pressure field.

Euler equations are typically represented in a "conserved" form, where the vector u contains (density, momentum, total energy) at each point x and time t. This function converts these conserved state variables to primitive variables (density, velocity, pressure) in order to calculate delta1 and delta2 according to [6].

The `u_values` parameter is a 3D array representing the conserved state variables at each point x and time t. This function iterates over each time step and calculates the primitive state vector for density, velocity, and pressure using appropriate transformations. The resulting primitive variables are stored in separate arrays `density_field`, `velocity_field`, and `pressure_field`, and returned.
"""
function convert_to_primitive(u_values)
    u_prim = zeros(size(u_values))
    for i in 1:size(u_values, 2)
        for j in 1:size(u_values, 3)
            # primitive_state_vector returns value without units
            u_p_M_T = primitive_state_vector(u_values[:, i, j]; gas=DRY_AIR)
            p_u = pressure(u_p_M_T[1], u_p_M_T[3]; gas=DRY_AIR)
            # Store density
            u_prim[1, i, j] = u_p_M_T[1]
            # Convert Mach to m/s using speed_of_sound
            u_prim[2, i, j] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[3]; gas=DRY_AIR)) 
            # Strip the unit of pressure so that it can be stored in an empty array
            u_prim[3, i, j] = ustrip(p_u)
        end
    end

    # Extract the density, velocity, and pressure fields
    density_field = u_prim[1, :, :]
    velocity_field = u_prim[2, :, :]
    pressure_field = u_prim[3, :, :]
    return density_field, velocity_field, pressure_field
end
