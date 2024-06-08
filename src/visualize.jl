using Plots
using FFMPEG
using ShockwaveProperties: primitive_state_vector, pressure, speed_of_sound, DRY_AIR
using Unitful: Pa, ustrip

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

# Euler is in a "conserved" form and the vector u contains (density, momentum, total energy) at each point x and time t
# So we want to convert these to (density, velocity, pressure) to calculate delta1 and delta2 according to [6]

# Read the data and create the animation
# x0_xmax, t_values, u_values, dims_u = read_output_file("C:/Users/user/Documents/School/Sem4/softwareentwicklungspraktikum/shock_wave_detection/ShockwaveProperties.jl/example/data/euler_scenario_2.out")

# Convert u_values to primitive variables. It's convert to density, velocity, pressure
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

function create_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)
    # Create a range of x values
    x = range(x0_xmax[1], stop=x0_xmax[2], length=dims_u[2])

    # Create a range of t values
    t = t_values

    # Create an animation
    anim = @animate for (i, t) in enumerate(t_values)
        p1 = Plots.plot(x, density_field[:, i], title="Density at Time $t", xlabel="x", ylabel="Density(kg/m^3)", label = "Density across x", size=(800, 600))
        p2 = Plots.plot(x, velocity_field[:, i], title="Velocity at Time $t", xlabel="x", ylabel="Velocity(m/s)", label = "Velocity across x", size=(800, 600))
        p3 = Plots.plot(x, pressure_field[:, i], title="Pressure at Time $t", xlabel="x", ylabel="Pressure(Pa)", label = "Pressure across x", size=(800, 600))
        plot(p1, p2, p3, layout = (3, 1))
    end
    # Save the animation as a gif
    gif(anim, "density_velocity_pressure_over_time.gif", fps = 10)
end

#density_field, velocity_field, pressure_field = convert_to_primitive(u_values)

#create_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)
