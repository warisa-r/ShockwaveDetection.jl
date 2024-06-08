include("visualize.jl")

# Use Interpolations.jl to perform np.gradient-like operations
using Interpolations

using StatsPlots

# Function to compute gradients using Interpolations.jl
# without explicit function definition
function compute_gradients(arr::AbstractVector, x)
    # Convert StepRangeLen to LinRange because that's what Interpolations.jl expects
    x_linrange = LinRange(first(x), last(x), length(x))
    itp = interpolate((x_linrange,), arr, Gridded(Linear()))
    grad = only.(Interpolations.gradient.(Ref(itp), x_linrange))
    return grad
end

# Function to compute zero crossings, detecting points where the second derivative changes sign\
# suggesting a potential shock
function zero_crossings(arr)
    return [i for i in 2:length(arr) if arr[i-1] * arr[i] < 0]
end

# Function to detect shock locations at a given time step
function detect_shocks_at_timestep(density_at_t, velocity_at_t, pressure_at_t, x)
    # Compute first gradients
    density_grad = compute_gradients(density_at_t, x)
    velocity_grad = compute_gradients(velocity_at_t, x)
    pressure_grad = compute_gradients(pressure_at_t, x)

    # Find points where the gradient exceeds a certain threshold
    threshold = 0.5  # Adjust this threshold as needed
    shock_location_density = findall(gradient -> abs(gradient) > threshold, density_grad)
    shock_location_velocity = findall(gradient -> abs(gradient) > threshold, velocity_grad)
    shock_location_pressure = findall(gradient -> abs(gradient) > threshold, pressure_grad)

    # Combine detections (common shock location across variables)
    shock_locations = intersect(intersect(shock_location_density, shock_location_velocity), shock_location_pressure)
    
    return shock_locations
end

# Read the data
x0_xmax, t_values, u_values, dims_u = read_output_file("C:/Users/user/Documents/School/Sem4/softwareentwicklungspraktikum/shock_wave_detection/ShockwaveProperties.jl/example/data/euler_scenario_2.out")
x = range(x0_xmax[1], stop=x0_xmax[2], length=dims_u[2])
density_field, velocity_field, pressure_field = convert_to_primitive(u_values)

# Detect shocks for each time step and store the positions
shock_positions_over_time = []

for t_step in 1: length(t_values)
    density_field_t = density_field[:, t_step]
    velocity_field_t = velocity_field[:, t_step]
    pressure_field_t = pressure_field[:, t_step]
    shock_positions = detect_shocks_at_timestep(density_field_t,velocity_field_t,pressure_field_t, x)
    push!(shock_positions_over_time, shock_positions)
end

# Create an animation
anim = @animate for (t_step, t) in enumerate(t_values)
    p1 = plot(x, density_field[:, t_step], title="Density at Time $t", xlabel="x", ylabel="Density", label = "Density across x", size=(800, 600))
    p2 = plot(x, velocity_field[:, t_step], title="Velocity at Time $t", xlabel="x", ylabel="Velocity", label = "Velocity across x", size=(800, 600))
    p3 = plot(x, pressure_field[:, t_step], title="Pressure at Time $t", xlabel="x", ylabel="Pressure", label = "Pressure across x", size=(800, 600))
    
    # Add markers for the shock positions
    shock_positions_t = shock_positions_over_time[t_step]
    for pos in shock_positions_t
        scatter!(p1, [x[pos]], [density_field[pos, t_step]], color=:red, label=false)
        scatter!(p2, [x[pos]], [velocity_field[pos, t_step]], color=:red, label=false)
        scatter!(p3, [x[pos]], [pressure_field[pos, t_step]], color=:red, label=false)
    end
    
    plot(p1, p2, p3, layout = (3, 1))
end

# Save the animation as a gif
gif(anim, "shock_over_time.gif", fps = 10)