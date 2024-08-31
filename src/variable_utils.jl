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

Euler equations are typically represented in a "conserved" form, where the vector u contains (density, momentum, total energy) at each point x and time t. This function converts these conserved state variables to primitive variables (density, velocity, pressure).

The `u_values` parameter is a 3D array representing the conserved state variables at each point x and time t. This function iterates over each time step and calculates the primitive state vector for density, velocity, and pressure using appropriate transformations. The resulting primitive variables are stored in separate arrays `density_field`, `velocity_field`, and `pressure_field`, and returned.
"""
function convert_to_primitive(u_values::Array{T, 3}, ncells, nsteps, mach_to_m_s=true) where T
    u_prim = similar(u_values)
    @threads for x in 1:ncells[1]
        for t in 1:nsteps
            # primitive_state_vector returns value without units
            u_p_M_T = primitive_state_vector(u_values[:, x, t], DRY_AIR)
            p_u = pressure(u_p_M_T[1], u_p_M_T[3], DRY_AIR)
            # Store density
            u_prim[1, x, t] = u_p_M_T[1]
            # Convert Mach to m/s using speed_of_sound
            if mach_to_m_s
                u_prim[2, x, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[3], DRY_AIR))
            else
                u_prim[2, x, t] = u_p_M_T[2]
            end
            # Strip the unit of pressure so that it can be stored in an empty array
            u_prim[3, x, t] = ustrip(p_u)
        end
    end

    # Extract the density, velocity, and pressure fields
    density_field = u_prim[1, :, :]
    velocity_field = u_prim[2, :, :]
    pressure_field = u_prim[3, :, :]
    return density_field, velocity_field, pressure_field
end

#TODO: make this more efficient. Either make it shorter or thread this??
function convert_to_primitive(u_values::Array{T, 4}, ncells, nsteps, mach_to_m_s=true) where T
    u_prim = similar(u_values)
    @threads for x in 1:ncells[1]
        for y in 1:ncells[2]
            for t in 1:nsteps
                # primitive_state_vector returns value without units
                u_p_M_T = primitive_state_vector(u_values[:, x, y, t], DRY_AIR)
                p_u = pressure(u_p_M_T[1], u_p_M_T[4], DRY_AIR)
                # Store density
                u_prim[1, x, y, t] = u_p_M_T[1]
                # Convert Mach to m/s using speed_of_sound
                if mach_to_m_s
                    u_prim[2, x, y, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR))
                    u_prim[3, x, y, t] = u_p_M_T[3] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR))
                else
                    u_prim[2, x, y, t] = u_p_M_T[2]
                    u_prim[3, x, y, t] = u_p_M_T[3]
                end
                # Strip the unit of pressure so that it can be stored in an empty array
                u_prim[4, x, y, t] = ustrip(p_u)
            end
        end
    end

    # Extract the density, velocity, and pressure fields
    density_field = u_prim[1, :, :, :]
    velocity_field = u_prim[2:3, :, :, :]
    pressure_field = u_prim[4, :, :, :]
    return density_field, velocity_field, pressure_field
end

"""
    cartesian_index_to_xy(shock_positions_t, x, y) -> Matrix

Converts Cartesian indices of shock positions into their corresponding x and y coordinates.

# Arguments
- `shock_positions_t`: An array of CartesianIndex indicating the positions of shocks.
- `x`: An array representing the x-coordinates in the domain.
- `y`: An array representing the y-coordinates in the domain.

# Returns
- `Matrix`: A 2xN matrix where the first row contains the x-coordinates and the second row contains the y-coordinates of the shock positions. The type of the elements in the matrix matches the type of elements in `x` and `y`.
"""
function cartesian_index_to_xy(shock_positions_t, x , y)
    coordinate_matrix = zeros(eltype(x), 2, length(shock_positions_t))
    # Extract x and y coordinates from CartesianIndex
    xs = [pos[1] for pos in shock_positions_t]
    ys = [pos[2] for pos in shock_positions_t]

    # Nail down where this is in x and y before scattering since shock_positions is just sets of indices not actual x and y values
    x_shocks = [x[x_pos] for x_pos in xs]
    y_shocks = [y[y_pos] for y_pos in ys]

    # Assign x_shocks to the first row and y_shocks to the second row of coordinate_matrix
    coordinate_matrix[1, :] = x_shocks
    coordinate_matrix[2, :] = y_shocks

    return coordinate_matrix
end