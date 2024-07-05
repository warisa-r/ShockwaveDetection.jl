# Since @animate is a macro from the Plots.jl package, we cannot specify importing it in the module.
using Plots
using CairoMakie: CairoMakie, heatmap!, record, Figure, Axis, vlines!, Colorbar, empty!
"""
    create_wave_animation(flow_data::FlowData)

Create an animation of density, velocity, and pressure fields over time from properties of the flow.

# Arguments
- `flow_data::FlowData`: A struct containing the flow data of x values, time steps, density field, velocity field, and pressure field.

# Returns
- This function saves an animation as a gif file named "density_velocity_pressure_over_time.gif" in the current directory and returns the animated object.
"""
function create_wave_animation(flow_data::FlowData; T=Float64)
    if typeof(flow_data.u) == Array{T, 3}
        # Handle the 1D flow case
        create_wave_animation_1D(flow_data::FlowData)
    elseif typeof(flow_data.u) == Array{T, 4}
        # Handle the 2D flow case
        create_wave_animation_2D(flow_data::FlowData)
    else
        error("Unsupported array dimensionality for flow_data.u")
    end
end

function create_wave_animation_1D(flow_data::FlowData)
    bounds = flow_data.bounds
    tsteps = flow_data.tsteps
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field

    # Create a range of x values
    x = range(bounds[1][1], stop=bounds[1][2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (i, t) in enumerate(tsteps)
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

function create_wave_animation_2D(flow_data::FlowData)
    #TODO: Implement this function
    # 6 Plots??? LOL would this even make any sense or is this even useful??
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
function create_wave_animation_with_shock(flow_data::FlowData, shock_positions_over_time; T=Float64)
    if typeof(flow_data.u) == Array{T, 3}
        # Handle the 1D flow case
        create_wave_animation_with_shock_1D(flow_data::FlowData, shock_positions_over_time)
    else
        error("Unsupported array dimensionality for flow_data.u")
    end
end

function create_wave_animation_with_shock_1D(flow_data::FlowData, shock_positions_over_time)
    bounds = flow_data.bounds
    tsteps = flow_data.tsteps
    density_field = flow_data.density_field
    velocity_field = flow_data.velocity_field
    pressure_field = flow_data.pressure_field

    # Create a range of x values
    x = range(bounds[1][1], stop=bounds[1][2], length=size(density_field, 1))

    # Create an animation
    anim = @animate for (t_step, t) in enumerate(tsteps)
        p1 = plot(x, density_field[:, t_step], title="Density at Time(kg/m^3) $t", xlabel="x", ylabel="Density", label="Density across x", size=(800, 600))
        if flow_data.mach_to_m_s
            p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time (m/s) $t", xlabel="x", ylabel="Velocity", label="Velocity across x", size=(800, 600))
        else
            p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time (Mach) $t", xlabel="x", ylabel="Velocity", label="Velocity across x", size=(800, 600))
        end
        p3 = plot(x, pressure_field[:, t_step], title="Pressure at Time(Pa) $t", xlabel="x", ylabel="Pressure", label="Pressure across x", size=(800, 600))
        
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

function create_heatmap_evo(flow_data::FlowData, field::Symbol; T = Float64)
    if typeof(flow_data.u) == Array{T, 3}
        # Handle the 1D flow case
        create_heatmap_evo_1D(flow_data, field)
    elseif typeof(flow_data.u) == Array{T, 4}
        # Handle the 2D flow case
        create_heatmap_evo_2D(flow_data, field)
    else
        error("Unsupported array dimensionality for flow_data.u")
    end
end

function create_heatmap_evo_1D(flow_data::FlowData, field::Symbol, tube_circumference=5.0)
    field_data = getfield(flow_data, field)
    tsteps = flow_data.tsteps
    bounds = flow_data.bounds

    # Define the y-range
    y = range(0.0, tube_circumference, length=7)

    # Define the x-range
    len_x = size(field_data, 1)
    x = range(bounds[1][1], bounds[1][2], length=len_x)

    # Create the figure for the animation
    fig = CairoMakie.Figure(size = (1000, 400))
    ax = CairoMakie.Axis(fig[1, 1], title = "$field Field Evolution")

    #TODO: add a colorbar

    # Record the animation
    CairoMakie.record(fig, "$(field)_evolution.gif", enumerate(tsteps); framerate = 10) do (t, t_step)
        field_t = field_data[:, t]
        field_tube = repeat(field_t, 1, 7)

        # Create the heatmap and store the returned object
        CairoMakie.heatmap!(ax, x, y, field_tube)
        
        ax.title = "$field Field - Time Step: $t_step"
    end
end

function create_heatmap_evo_2D(flow_data, field)
    if field == :velocity_field
        velocity_field = getfield(flow_data, field)
        velocity_field_x = velocity_field[1, :, :, :]
        velocity_field_y = velocity_field[2, :, :, :]
        field_data = sqrt.(velocity_field_x.^2 + velocity_field_y.^2)
    else
        field_data = getfield(flow_data, field)
    end

    tsteps = flow_data.tsteps
    bounds = flow_data.bounds
    ncells = flow_data.ncells

    # Define the x-range
    x = range(bounds[1][1], bounds[1][2], length=ncells[1])

    # Define the y-range
    y = range(bounds[2][1], bounds[2][2], length=ncells[2])

    # Create the figure for the animation
    fig = CairoMakie.Figure(size = (1000, 800))
    ax = CairoMakie.Axis(fig[1, 1], title = "$field Field Evolution")

    # Record the animation
    CairoMakie.record(fig, "$(field)_evolution.gif", enumerate(tsteps); framerate = 10) do (t, t_step)
        
        # Extract the field data for the current time step
        field_t = field_data[:, :, t]

        # Create the heatmap and store the returned object
        CairoMakie.heatmap!(ax, x, y, field_t)
    
        # Update the title with the current time step
        ax.title = "$field Field - Time Step: $t_step"
    end
end

function create_heatmap_evo_with_shock(flow_data, detection, field::Symbol = :density_field, show_normal_vector = true; T= Float64)
    if typeof(flow_data.u) == Array{T, 3}
        # Handle the 1D flow case
        #TODO: pipeline of detection in 1D case hasnt been implemented yet. It should contain the field shock_positions_over_time as well!
        shock_positions_over_time = detection.shock_positions_over_time
        create_heatmap_evo_with_shock_1D(flow_data, shock_positions_over_time, field)
    elseif typeof(flow_data.u) == Array{T, 4}
        # Handle the 2D flow case
        shock_positions_over_time = detection.shock_positions_over_time
        shock_fits_over_time = detection.shock_fits_over_time
        create_heatmap_evo_with_shock_2D(flow_data, detection, field, show_normal_vector)
    else
        error("Unsupported array dimensionality for flow_data.u")
    end
end

function create_heatmap_evo_with_shock_1D(flow_data::FlowData, shock_positions_over_time, field::Symbol, tube_circumference=5.0)
    field_data = getfield(flow_data, field)
    tsteps = flow_data.tsteps
    bounds = flow_data.bounds

    # Define the "fake" y-range
    y = range(0.0, tube_circumference, length=7)

    # Define the x-range
    len_x = size(field_data, 1)
    x = range(bounds[1][1], bounds[1][2], length=len_x)

    # Create the figure for the animation
    fig = CairoMakie.Figure(size = (1000, 400))
    ax = CairoMakie.Axis(fig[1, 1], title = "$field Field Evolution")

    # Record the animation
    CairoMakie.record(fig, "$(field)_evolution.gif", enumerate(tsteps); framerate = 10) do (t, t_step)
        
        # Extract the field data for the current time step
        field_t = field_data[:, t]
        field_tube = repeat(field_t, 1, 7)

        # Create the heatmap and store the returned object
        CairoMakie.heatmap!(ax, x, y, field_tube)
    
        # Update the title with the current time step
        ax.title = "$field Field - Time Step: $t_step"
        
        # Plot the positions of the shock waves
        shock_positions = shock_positions_over_time[t]
        x_shocks = [x[pos] for pos in shock_positions]
        if length(shock_positions) > 0
            CairoMakie.vlines!(ax, x_shocks, color = :red, linestyle = :solid, linewidth = 2.0)
        end
    end
end

function create_heatmap_evo_with_shock_2D(flow_data, detection, field, show_normal_vector)
    if field == :velocity_field
        velocity_field = getfield(flow_data, field)
        velocity_field_x = velocity_field[1, :, :, :]
        velocity_field_y = velocity_field[2, :, :, :]
        field_data = sqrt.(velocity_field_x.^2 + velocity_field_y.^2)
    else
        field_data = getfield(flow_data, field)
    end

    
    tsteps = flow_data.tsteps
    bounds = flow_data.bounds
    ncells = flow_data.ncells

    shock_clusters_over_time = detection.shock_clusters_over_time
    shock_positions_over_time = detection.shock_positions_over_time
    shock_fits_over_time = detection.shock_fits_over_time


    # Define the x-range
    x = range(bounds[1][1], bounds[1][2], length=ncells[1])

    # Define the y-range
    y = range(bounds[2][1], bounds[2][2], length=ncells[2])

    # Create the figure for the animation
    fig = CairoMakie.Figure(size = (1000, 800))
    ax = CairoMakie.Axis(fig[1, 1], title = "$field Field Evolution")

    # Set explicit limits for the plot axes
    CairoMakie.xlims!(ax, bounds[1][1], bounds[1][2])
    CairoMakie.ylims!(ax, bounds[2][1], bounds[2][2])

    # Record the animation
    CairoMakie.record(fig, "$(field)_evolution.gif", enumerate(tsteps); framerate = 10) do (t, t_step)
        shock_clusters = shock_clusters_over_time[t]
        num_clusters = length(shock_clusters)
        
        # Extract the field data for the current time step
        field_t = field_data[:, :, t]

        # Create the heatmap and store the returned object
        CairoMakie.heatmap!(ax, x, y, field_t)

        # Overlay shock positions if available for the current time step
        if !isempty(shock_positions_over_time[t])
            shock_positions = shock_positions_over_time[t]
            # Extract x and y coordinates from CartesianIndex
            xs = [pos[1] for pos in shock_positions]
            ys = [pos[2] for pos in shock_positions]

            # Nail down where this is in x and y before scattering since shock_positions is just sets of indices not actual x and y values
            x_shocks = [x[x_pos] for x_pos in xs]
            y_shocks = [y[y_pos] for y_pos in ys]
            CairoMakie.scatter!(ax, x_shocks, y_shocks, color = :red, markersize = 3)

            show_curve = true

            if show_curve
                shock_fits = shock_fits_over_time[t]
                for shock_fit in shock_fits
                    if shock_fit.model == vline_model
                        CairoMakie.vlines!(ax, shock_fit.parameters[1], color=:orange)
                        if show_normal_vector
                            average_y = bounds[2][1] + bounds[2][2]/ 2
                            start_x = shock_fit.parameters[1]
                            # Here evenly_spaced_range vector can be an empty set because one normal vector is enough to represent
                            # the vertical line
                            normals_x , normals_y = calculate_normal_vector(shock_fit, [], flow_data, t)
                            CairoMakie.arrows!(ax, [start_x], [average_y], normals_x, normals_y, color=:green)
                        end
                    else
                        if shock_fit.model == circle_model
                            angles = range(shock_fit.range[1], shock_fit.range[2], length= round(Int, ncells[1] / num_clusters))
                            # Calculate x and y coordinates based on the circle equation
                            x_values = shock_fit.parameters[1] .+ shock_fit.parameters[3] .* cos.(angles)
                            y_values = shock_fit.parameters[2] .+ shock_fit.parameters[3] .* sin.(angles)
                            if show_normal_vector
                                normals_x, normals_y = calculate_normal_vector(shock_fit, angles, flow_data, t)
                            end
                            
                        else
                            x_values = range(shock_fit.range[1], shock_fit.range[2], length= round(Int, ncells[1] / num_clusters))
                            if shock_fit.model == line_model
                                y_values = shock_fit.parameters[1] .+ shock_fit.parameters[2] .* x_values
                            elseif shock_fit.model == parabola_model
                                y_values = shock_fit.parameters[1] .* x_values.^2 .+ shock_fit.parameters[2] .* x_values .+ shock_fit.parameters[3]
                            end
                            if show_normal_vector
                                normals_x, normals_y = calculate_normal_vector(shock_fit, x_values, flow_data, t)
                            end
                        end
                        if show_normal_vector
                            CairoMakie.arrows!(ax, x_values, y_values, normals_x, normals_y, color=:green) # Normal vector of the circle
                        end
                        CairoMakie.lines!(ax, x_values, y_values, color=:orange)
                    end
                end
            end
        end

        println("Frame rendered for time step: $t_step")
    
        # Update the title with the current time step
        ax.title = "$field Field - Time Step: $t_step"
    end
end

function plot_shock_fits(flow_data, shock_clusters, fits, show_normal_vector, t)
    bounds = flow_data.bounds
    ncells = flow_data.ncells

    num_clusters = length(shock_clusters)

    # Create the figure for the animation
    fig = CairoMakie.Figure(size = (1000, 800))
    ax = CairoMakie.Axis(fig[1, 1], title = "Plot of the fitted curve")

    # Set explicit limits for the plot axes
    CairoMakie.xlims!(ax, bounds[1][1], bounds[1][2])
    CairoMakie.ylims!(ax, bounds[2][1], bounds[2][2])

    # If the cluster is empty then there is also no fit to plot
    if isempty(shock_clusters)
        return fig
    end

    for shock_cluster in shock_clusters
        if !isempty(shock_cluster)
            # Extract x and y coordinates
            x_shocks = [point[1] for point in shock_cluster]
            y_shocks = [point[2] for point in shock_cluster]

            CairoMakie.scatter!(ax, x_shocks, y_shocks, color = :red, markersize = 3)
        end
    end

    if isempty(fits)
        return fig
    end

    for fit in fits
        if fit.model == vline_model
            CairoMakie.vlines!(ax, fit.parameters[1], color=:blue)
            if show_normal_vector
                # Calculate normal vectors
                # Plot the normal vector as an arrow at the average of y-values for visibility
                average_y = bounds[2][1] + bounds[2][2]/ 2
                start_x = fit.parameters[1]
                # Here evenly_spaced_range vector can be an empty set because one normal vector is enough to represent
                # the vertical line
                normals_x , normals_y = calculate_normal_vector(fit, [], flow_data, t)
                CairoMakie.arrows!(ax, [start_x], [average_y], normals_x, normals_y, color=:green)
            end
        else
            if fit.model == circle_model
                angles = range(fit.range[1], fit.range[2], length= round(Int, ncells[1] / num_clusters))
                # Calculate x and y coordinates based on the circle equation
                x_values = fit.parameters[1] .+ fit.parameters[3] .* cos.(angles)
                y_values = fit.parameters[2] .+ fit.parameters[3] .* sin.(angles)
                if show_normal_vector
                    normals_x, normals_y = calculate_normal_vector(fit, angles, flow_data, t)
                end
                
            else
                x_values = range(fit.range[1], fit.range[2], length= round(Int, ncells[1] / num_clusters))
                if fit.model == line_model
                    y_values = fit.parameters[1] .+ fit.parameters[2] .* x_values
                elseif fit.model == parabola_model
                    y_values = fit.parameters[1] .* x_values.^2 .+ fit.parameters[2] .* x_values .+ fit.parameters[3]
                end
                if show_normal_vector
                    normals_x, normals_y = calculate_normal_vector(fit, x_values, flow_data, t)
                end
            end
            if show_normal_vector
                CairoMakie.arrows!(ax, x_values, y_values, normals_x, normals_y, color=:green) # Normal vector of the circle
            end
            CairoMakie.lines!(ax, x_values, y_values, color=:blue)
        end
    end

    return fig
end

"""
    plot_shock_fits_over_time(flow_data, detection, show_normal_vector = false; T=Float64)

Plot the shock fits over time for a given flow data and detection results. This differs from other visualization function in that it shows the shock fits and doesn't save it as gif file.

# Arguments
- `flow_data`: The flow data containing information about the flow field.
- `detection`: An object containing the shock detection results, including shock clusters and fits over time.
- `show_normal_vector`: A boolean flag indicating whether to show the normal vector for each shock fit. Defaults to `false`.
- `T`: The type of the elements in the flow data's velocity field. Defaults to `Float64`.

"""
function plot_shock_fits_over_time(flow_data, detection, show_normal_vector = false; T=Float64)
    nsteps = flow_data.nsteps

    shock_clusters_over_time = detection.shock_clusters_over_time
    shock_fits_over_time = detection.shock_fits_over_time
    
    if typeof(flow_data.u) == Array{T, 3}
        println("Feature doesn't support 1D case")
    else
        for t in 1:nsteps
            fig = plot_shock_fits(flow_data, shock_clusters_over_time[t], shock_fits_over_time[t], show_normal_vector, t)
            display(fig)
        end
    end
    
end