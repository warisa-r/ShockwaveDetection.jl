# Since @animate is a macro from the Plots.jl package, we cannot specify importing it in the module.
using Plots

"""
    create_wave_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)

Create an animation of density, velocity, and pressure fields over time.

# Arguments
- `x0_xmax::Tuple{Float64, Float64}`: A tuple containing the minimum and maximum x values.
- `t_values::Array{Float64}`: An array of time values.
- `density_field::Array{Float64, 2}`: A 2D array representing the density field. Each column represents the density at a different time.
- `velocity_field::Array{Float64, 2}`: A 2D array representing the velocity field. Each column represents the velocity at a different time.
- `pressure_field::Array{Float64, 2}`: A 2D array representing the pressure field. Each column represents the pressure at a different time.

# Returns
- This function saves an animation as a gif file named "density_velocity_pressure_over_time.gif" in the current directory and returns nothing.
"""
function create_wave_animation(x0_xmax, t_values, density_field, velocity_field, pressure_field)
    # Create a range of x values
    x = range(x0_xmax[1], stop=x0_xmax[2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (i, t) in enumerate(t_values)
        p1 = plot(x, density_field[:, i], title="Density at Time $t", xlabel="x", ylabel="Density(kg/m^3)", label="Density across x", size=(800, 600))
        p2 = plot(x, velocity_field[:, i], title="Velocity at Time $t", xlabel="x", ylabel="Velocity(m/s)", label="Velocity across x", size=(800, 600))
        p3 = plot(x, pressure_field[:, i], title="Pressure at Time $t", xlabel="x", ylabel="Pressure(Pa)", label="Pressure across x", size=(800, 600))
        plot(p1, p2, p3, layout=(3, 1))
    end

    # Save the animation as a gif
    gif(anim, "density_velocity_pressure_over_time.gif", fps=10)

    return anim
end

"""
    create_wave_animation_with_shock(x, t_values, density_field, velocity_field, pressure_field, shock_positions_over_time; save_file=false)

Create an animation of density, velocity, and pressure fields over time, with markers for shock positions.

# Arguments
- `x0_xmax::Tuple{Float64, Float64}`: A tuple containing the minimum and maximum x values.
- `t_values::Vector`: Time values.
- `density_field::Matrix`: Density field over space and time.
- `velocity_field::Matrix`: Velocity field over space and time.
- `pressure_field::Matrix`: Pressure field over space and time.
- `shock_positions_over_time::Vector{Vector{Int}}`: Shock positions over time.
- `save_file::Bool`: Whether to save the animation as a gif file (default: false).

# Returns
- `anim::Animation`: Animation object.

This function creates an animation of density, velocity, and pressure fields over time, with markers indicating shock positions at each time step. The animation is generated using the Plots.jl package and returned as an Animation object. If `save_file` is true, the animation is saved as a gif file with the filename 'density_velocity_pressure_over_time_with_shock_positions.gif'.
"""
function create_wave_animation_with_shock(x0_xmax, t_values, density_field, velocity_field, pressure_field, shock_positions_over_time)
    # Create a range of x values
    x = range(x0_xmax[1], stop=x0_xmax[2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (t_step, t) in enumerate(t_values)
        p1 = plot(x, density_field[:, t_step], title="Density at Time $t", xlabel="x", ylabel="Density", label="Density across x", size=(800, 600))
        p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time $t", xlabel="x", ylabel="Velocity", label="Velocity across x", size=(800, 600))
        p3 = plot(x, pressure_field[:, t_step], title="Pressure at Time $t", xlabel="x", ylabel="Pressure", label="Pressure across x", size=(800, 600))
        
        # Add markers for the shock positions
        shock_positions_t = shock_positions_over_time[t_step]
        for pos in shock_positions_t
            scatter!(p1, [x[pos]], [density_field[pos, t_step]], color=:red, label=false)
            scatter!(p2, [x[pos]], [velocity_field[pos, t_step]], color=:red, label=false)
            scatter!(p3, [x[pos]], [pressure_field[pos, t_step]], color=:red, label=false)
        end
        
        plot(p1, p2, p3, layout=(3, 1))
    end

    gif(anim, "density_velocity_pressure_over_time_with_shock_positions.gif", fps=10)

    return anim
end