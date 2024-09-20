using LsqFit
using Plots  # Import the Plots library
using BenchmarkTools  # Import the BenchmarkTools library
using Interpolations
using LinearAlgebra

# circle model function
function circle_model(xy, p)
    x0, y0, r = p
    x = xy[:, 1]
    y = xy[:, 2]
    return (x .- x0).^2 + (y .- y0).^2 .- r^2
end

function fit_circle_with_random_guess(xy)
    # Random initial guess for circle fitting
    p0 = rand(3)  # Random initial guesses for a, b, c
    println("Random initial guess: ", p0)

    # Perform the fit
    fit = curve_fit(circle_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Improved function to fit a circle with a better initial guess
function fit_circle_with_improved_guess(xy)
    
    # Improved initial guess for circle fitting
    x = xy[:, 1]
    y = xy[:, 2]

    p0 = [mean(x), mean(y), sqrt(mean((x .- mean(x)).^2 .+ (y .- mean(y)).^2))]    

    # Perform the fit
    fit = curve_fit(circle_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Function to generate noisy circular data
function generate_noisy_circular_data(p, noise_level)
    x0, y0, r = p
    num_points = 100

    θ = LinRange(0, 2π, num_points)
    x = x0 .+ r .* cos.(θ) .+ noise_level .* randn(num_points)
    y = y0 .+ r .* sin.(θ) .+ noise_level .* randn(num_points)
    # Add noise
    y_noisy = y .+ randn(length(y)) * noise_level  # Adjust noise level as needed

    # Combine into xy matrix
    xy = hcat(x, y_noisy)
    return xy
end

# Define the parameter combinations to test
# Define the parameter combinations to test
param_combinations = [
    [49, -49, 0.1],
    [0.005, 0.005, 10],
    [1, 1, 20],
    [25, 25, 15],
    [-25, -25, 30],
    [50, 50, 5],
    [-50, -50, 10]
]

noise_levels = [0.0, 1.0, 5.0, 7.5, 10.0]

for (i, params) in enumerate(param_combinations)
    results = []
    for noise_level in noise_levels
        xy = generate_noisy_circular_data(params, noise_level)
        
        # Benchmark improved initial guess
        improved_benchmark = @benchmark fit_circle_with_improved_guess($xy)
        improved_initial_guess, improved_fit, _ = fit_circle_with_improved_guess(xy)
        improved_error = norm(improved_fit.param - params)
        
        # Benchmark random initial guess
        random_benchmark = @benchmark fit_circle_with_random_guess($xy)
        random_initial_guess, random_fit, _ = fit_circle_with_random_guess(xy)
        random_error = norm(random_fit.param - params)
        
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
    savefig(p, "sandbox/param_performance_results/circle/benchmark_results_params_$(i)_noise.png")
end