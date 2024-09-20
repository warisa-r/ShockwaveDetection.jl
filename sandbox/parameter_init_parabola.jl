using LsqFit
using Plots  # Import the Plots library
using BenchmarkTools  # Import the BenchmarkTools library
using Interpolations
using LinearAlgebra

# Parabola model function
function parabola_model(xy, p)
    a, b, c = p
    x = xy[:, 1]
    y = xy[:, 2]
    return a * x.^2 .+ b * x .+ c .- y
end

function fit_parabola_with_random_guess(xy)
    # Random initial guess for parabola fitting
    p0 = rand(3)  # Random initial guesses for a, b, c
    println("Random initial guess: ", p0)

    # Perform the fit
    fit = curve_fit(parabola_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Improved function to fit a parabola with a better initial guess
function fit_parabola_with_improved_guess(xy)
    
    # Improved initial guess for parabola fitting
    x_range = LinRange(-50, 50, 100)
    y_range = LinRange(-50, 50, 100)
    x = xy[:, 1]
    y = xy[:, 2]

    c_guess = mean(y)

    itp = interpolate((y_range,), y, Gridded(Linear()))
    # Evaluate the gradient at the points in xy
    first_derivative = only.(gradient.(Ref(itp), y_range))

    b_guess = mean(first_derivative)
    println("b approximate: ", b_guess)

    # Create an interpolation object for the first derivative
    itp_first_derivative = interpolate((y_range,), first_derivative, Gridded(Linear()))

    # Evaluate the second derivative at the points in x
    second_derivative = only.(gradient.(Ref(itp_first_derivative), y_range))
    # Compute the mean of the second derivative
    mean_second_derivative = mean(second_derivative)
    a_guess = mean_second_derivative/2
    println("a approximate (interpolation guess): ", a_guess)

    p0 = [a_guess, b_guess, c_guess]

    # Perform the fit
    fit = curve_fit(parabola_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Function to generate noisy parabolic data
function generate_noisy_parabolic_data(p, noise_level)
    a, b, c = p
    # Generate data points along a parabola
    x = range(-50, 50, length=100)
    y = a .* x.^2 + b .* x .+ c

    # Add noise
    y_noisy = y .+ randn(length(y)) * noise_level  # Adjust noise level as needed

    # Combine into xy matrix
    xy = hcat(x, y_noisy)
    return xy
end

# Define the parameter combinations to test
param_combinations = [
    [0.005, 0.005, 0.005],
    [1, 1, 1],
    [500, 500, 500],
    [0.005, 500, 1],
    [500, 0.005, 1],
    [0.005, 5000, 0.5]
]

noise_levels = [0.0, 0.1, 1.0, 10.0]

for (i, params) in enumerate(param_combinations)
    results = []
    for noise_level in noise_levels
        xy = generate_noisy_parabolic_data(params, noise_level)
        
        # Benchmark improved initial guess
        improved_benchmark = @benchmark fit_parabola_with_improved_guess($xy)
        improved_initial_guess, improved_fit, improved_error = fit_parabola_with_improved_guess(xy)
        
        # Benchmark random initial guess
        random_benchmark = @benchmark fit_parabola_with_random_guess($xy)
        random_initial_guess, random_fit, random_error = fit_parabola_with_random_guess(xy)
        
        push!(results, (noise_level, improved_benchmark, improved_error, random_benchmark, random_error))
    end

    # Extract time and memory allocation data
    improved_times = [minimum(r[2]).time / 1e9 for r in results] # Convert to seconds
    improved_memories = [minimum(r[2]).memory / 1024 for r in results] # Convert to KB
    random_times = [minimum(r[4]).time / 1e9 for r in results] # Convert to seconds
    random_memories = [minimum(r[4]).memory / 1024 for r in results] # Convert to KB
    improved_errors = [r[3] for r in results]
    random_errors = [r[5] for r in results]

    # Plot the results once for all noise levels
    p = plot(layout = (3, 1), size = (800, 1200))
    
    # Time plot
    plot!(p[1], noise_levels, improved_times, label="Improved Time (s)", xlabel="Noise Level", ylabel="Time (s)", legend=:topright, title="Time Comparison (Params: $params)")
    plot!(p[1], noise_levels, random_times, label="Random Time (s)", xlabel="Noise Level", ylabel="Time (s)", legend=:topright)
    
    # Memory plot
    plot!(p[2], noise_levels, improved_memories, label="Improved Memory (KB)", xlabel="Noise Level", ylabel="Memory (KB)", legend=:topright, title="Memory Comparison (Params: $params)")
    plot!(p[2], noise_levels, random_memories, label="Random Memory (KB)", xlabel="Noise Level", ylabel="Memory (KB)", legend=:topright)
    
    # Error plot
    plot!(p[3], noise_levels, improved_errors, label="Improved Error", xlabel="Noise Level", ylabel="Error", legend=:topright, title="Error Comparison (Params: $params)")
    plot!(p[3], noise_levels, random_errors, label="Random Error", xlabel="Noise Level", ylabel="Error", legend=:topright)

    # Save the plot to a file
    savefig(p, "sandbox/param_performance_results/parabola/benchmark_results_params_$(i)_noise.png")
end