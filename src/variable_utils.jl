using ShockwaveProperties: primitive_state_vector, pressure, speed_of_sound, DRY_AIR
using Euler2D: Euler2D
using Unitful: ustrip
using Distributions

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
function convert_to_primitive(u_values::Array{T, 3}, ncells, nsteps, mach_to_m_s=true, noise=true) where T
    u_prim = similar(u_values)
    noise_dist = Normal(0.0, 0.01)  # Gaussian noise with mean 0 and standard deviation as noise_level
    @threads for x in 1:ncells[1]
        for t in 1:nsteps
            # primitive_state_vector returns value without units
            u_p_M_T = primitive_state_vector(u_values[:, x, t], DRY_AIR)
            p_u = pressure(u_p_M_T[1], u_p_M_T[3], DRY_AIR)
            # Store density with/without added Gaussian noise
            if noise
                u_prim[1, x, t] = u_p_M_T[1] + rand(noise_dist)
            else
                u_prim[1, x, t] = u_p_M_T[1]
            end
            # Convert Mach to m/s using speed_of_sound and add noise
            if mach_to_m_s
                if noise
                    u_prim[2, x, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[3], DRY_AIR)) + rand(noise_dist)
                else 
                    u_prim[2, x, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[3], DRY_AIR))
                end
            else
                if noise
                    u_prim[2, x, t] = u_p_M_T[2] + rand(noise_dist)
                else 
                    u_prim[2, x, t] = u_p_M_T[2]
                end
            end
            # Strip the unit of pressure so that it can be stored in an empty array and add noise
            u_prim[3, x, t] = ustrip(p_u) + rand(noise_dist)
        end
    end

    # Extract the density, velocity, and pressure fields
    density_field = u_prim[1, :, :]
    velocity_field = u_prim[2, :, :]
    pressure_field = u_prim[3, :, :]
    return density_field, velocity_field, pressure_field
end

function convert_to_primitive(u_values::Array{T, 4}, ncells, nsteps, mach_to_m_s=true, noise=true) where T
    u_prim = similar(u_values)
    noise_dist = Normal(0.0, 0.1)
    @threads for x in 1:ncells[1]
        for y in 1:ncells[2]
            for t in 1:nsteps
                # primitive_state_vector returns value without units
                u_p_M_T = primitive_state_vector(u_values[:, x, y, t], DRY_AIR)
                p_u = pressure(u_p_M_T[1], u_p_M_T[4], DRY_AIR)
                # Store density with/without noise
                if noise
                    u_prim[1, x, y, t] = u_p_M_T[1] + rand(noise_dist)
                else 
                    u_prim[1, x, y, t] = u_p_M_T[1]
                end
                # Convert Mach to m/s using speed_of_sound and add noise
                if mach_to_m_s
                    if noise
                        u_prim[2, x, y, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR)) + rand(noise_dist)
                        u_prim[3, x, y, t] = u_p_M_T[3] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR)) + rand(noise_dist)
                    else 
                        u_prim[2, x, y, t] = u_p_M_T[2] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR)) 
                        u_prim[3, x, y, t] = u_p_M_T[3] * ustrip(speed_of_sound(u_p_M_T[4], DRY_AIR)) 
                    end
                else
                    if noise
                        u_prim[2, x, y, t] = u_p_M_T[2] + rand(noise_dist)
                        u_prim[3, x, y, t] = u_p_M_T[3] + rand(noise_dist)
                    else 
                        u_prim[2, x, y, t] = u_p_M_T[2] 
                        u_prim[3, x, y, t] = u_p_M_T[3] 
                    end
                end
                # Strip the unit of pressure so that it can be stored in an empty array (and add noise)
                if noise
                    u_prim[4, x, y, t] = ustrip(p_u) + rand(noise_dist)
                else 
                    u_prim[4, x, y, t] = ustrip(p_u)
                end
            end
        end
    end

    # Extract the density, velocity, and pressure fields
    density_field = u_prim[1, :, :, :]
    velocity_field = u_prim[2:3, :, :, :]
    pressure_field = u_prim[4, :, :, :]
    return density_field, velocity_field, pressure_field
end

function convert_to_primitive(sim_data, nsteps, mach_to_m_s=false, noise=true)
    density_field = []
    velocity_field = []
    pressure_field = []
    
    noise_dist = Normal(0.0, 0.1)  # Gaussian noise with mean 0 and standard deviation as noise_level

    @threads for t in 1:nsteps
        if noise 
            density_t = [x !== nothing ? ustrip(x) + rand(noise_dist) : NaN for x in Euler2D.density_field(sim_data, t)]
            pressure_t = Euler2D.pressure_field(sim_data, t, DRY_AIR)
            velocity_t = [x !== nothing ? ustrip(x) + rand(noise_dist) : NaN for x in Euler2D.velocity_field(sim_data, t)]
        else 
            density_t = [x !== nothing ? ustrip(x) : NaN for x in Euler2D.density_field(sim_data, t)]
            pressure_t = Euler2D.pressure_field(sim_data, t, DRY_AIR)
            velocity_t = [x !== nothing ? ustrip(x) : NaN for x in Euler2D.velocity_field(sim_data, t)]
        end
       # TODO: find a way to make this work properly. Since speed of sound needs density or pressure and we have no access to temperature 
        if mach_to_m_s
            if noise 
                speed_of_sound_t = [x !== nothing ? ustrip(speed_of_sound(ustrip(x), DRY_AIR)) + rand(noise_dist) : NaN for x in pressure_t]
                #velocity_t = velocity_t .* speed_of_sound_t
            else 
                speed_of_sound_t = [x !== nothing ? ustrip(speed_of_sound(ustrip(x), DRY_AIR)) : NaN for x in pressure_t]
            end
        end
        if noise
            pressure_t = [x !== nothing ? ustrip(x) + rand(noise_dist) : NaN for x in pressure_t]
        else
            pressure_t = [x !== nothing ? ustrip(x) : NaN for x in pressure_t]
        end

        push!(density_field, density_t)
        push!(pressure_field, pressure_t)
        push!(velocity_field, velocity_t)
    end

    # Convert lists of lists into multidimensional arrays
    density_field_array = cat(density_field..., dims=3)
    velocity_field_array = cat(velocity_field..., dims=4)
    pressure_field_array = cat(pressure_field..., dims=3)

    return density_field_array, velocity_field_array, pressure_field_array
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