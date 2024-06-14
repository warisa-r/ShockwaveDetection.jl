# Since @animate is a macro from the Plots.jl package, we cannot specify importing it in the module.
using Plots

"""
    create_wave_animation(flow_data::FlowData)

Create an animation of density, velocity, and pressure fields over time from properties of the flow.

# Arguments
- `flow_data::FlowData`: A struct containing the flow data of x values, time steps, density field, velocity field, and pressure field.

# Returns
- This function saves an animation as a gif file named "density_velocity_pressure_over_time.gif" in the current directory and returns the animated object.
"""
function create_wave_animation(flow_data::FlowData)
    x0_xmax = flow_data.x0_xmax
    t_values = flow_data.t_values
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field

    # Create a range of x values
    x = range(x0_xmax[1], stop=x0_xmax[2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (i, t) in enumerate(t_values)
        p1 = plot(x, density_field[:, i], title="Density at Time $t", xlabel="x", ylabel="Density(kg/m^3)", label="Density across x", size=(800, 600))
        if flow_data.mach_to_m_s
            p2 = plot(x, velocity_field[:, i], title="Velocity at Time $t", xlabel="x", ylabel="Velocity(m/s)", label="Velocity across x", size=(800, 600))
        else
            p2 = plot(x, velocity_field[:, i], title="Velocity at Time $t", xlabel="x", ylabel="Velocity(Mach)", label="Velocity across x", size=(800, 600))
        end
        
        p3 = plot(x, pressure_field[:, i], title="Pressure at Time $t", xlabel="x", ylabel="Pressure(Pa)", label="Pressure across x", size=(800, 600))
        plot(p1, p2, p3, layout=(3, 1))
    end

    # Save the animation as a gif
    gif(anim, "density_velocity_pressure_over_time.gif", fps=10)

    return anim
end

"""
    create_wave_animation_with_shock(flow_data::FlowData, shock_positions_over_time)

Create an animation of density, velocity, and pressure fields over time from properties of the flow with shock wave positions.

# Arguments
- `flow_data::FlowData`: A struct containing the flow data of x values, time steps, density field, velocity field, and pressure field.
- `shock_positions_over_time::Vector{Vector{Int}}`: A vector of vectors containing the shock positions over time from the detection algorithm.

# Returns
- This function saves an animation as a gif file named "density_velocity_pressure_over_time_with_shock_positions.gif" in the current directory 
and returns the animated object.
"""
function create_wave_animation_with_shock(flow_data::FlowData, shock_positions_over_time)
    x0_xmax = flow_data.x0_xmax
    t_values = flow_data.t_values
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field

    # Create a range of x values
    x = range(x0_xmax[1], stop=x0_xmax[2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (t_step, t) in enumerate(t_values)
        p1 = plot(x, density_field[:, t_step], title="Density at Time(kg/m^3) $t", xlabel="x", ylabel="Density", label="Density across x", size=(800, 600))
        if flow_data.mach_to_m_s
            p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time (m/s) $t", xlabel="x", ylabel="Velocity", label="Velocity across x", size=(800, 600))
        else
            p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time (Mach) $t", xlabel="x", ylabel="Velocity", label="Velocity across x", size=(800, 600))
        end
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