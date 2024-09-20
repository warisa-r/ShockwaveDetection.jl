using LsqFit
using Plots  # Import the Plots library
using BenchmarkTools  # Import the BenchmarkTools library
using Interpolations
using LinearAlgebra

# Horizontal line model function (only intercept b)
function hline_model(xy, p)
    b = p[1]
    y = xy[:, 2]
    return y .- b
end

# Function to fit a horizontal line with a random initial guess
function fit_hline_with_random_guess(xy)
    p0 = rand(1)  # Random initial guess for b
    println("Random initial guess: ", p0)

    fit = curve_fit(hline_model, xy, zeros(size(xy, 1)), p0)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Improved function to fit a horizontal line with a better initial guess
function fit_hline_with_improved_guess(xy)
    y = xy[:, 2]

    b_guess = mean(y)  # The best guess for horizontal line b is the mean of y
    p0 = [b_guess]

    fit = curve_fit(hline_model, xy, zeros(size(xy, 1)), p0)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Function to generate noisy horizontal data
function generate_noisy_hline_data(p, noise_level)
    b = p[1]
    x = range(-50, 50, length=100)
    y = b .+ randn(length(x)) * noise_level  # Just noise around the intercept b

    xy = hcat(x, y)
    return xy
end

# Define the parameter combinations to test (just the intercept b)
param_combinations = [
    [0.5],
    [1.0],
    [500.0],
    [0.005],
    [5000.0],
    [0.0005]
]

noise_levels = [0.0, 1.0, 5.0, 7.5, 10.0]

# Define the output path
output_path = raw"C:\Users\eivan\Dropbox\My PC (DESKTOP-AB2UFDR)\Desktop\SEP_Shockwave\1\param_performance_results\hline"

# Create the directory if it doesn't exist
isdir(output_path) || mkpath(output_path)

for (i, params) in enumerate(param_combinations)
    results = []
    for noise_level in noise_levels
        xy = generate_noisy_hline_data(params, noise_level)

        improved_benchmark = @benchmark fit_hline_with_improved_guess($xy)
        improved_initial_guess, improved_fit, _ = fit_hline_with_improved_guess(xy)
        improved_error = norm(improved_fit.param - params)

        random_benchmark = @benchmark fit_hline_with_random_guess($xy)
        random_initial_guess, random_fit, _ = fit_hline_with_random_guess(xy)
        random_error = norm(random_fit.param - params)

        push!(results, (noise_level, improved_benchmark, improved_error, random_benchmark, random_error))
    end

    improved_times = [minimum(r[2]).time / 1e9 for r in results]
    improved_memories = [minimum(r[2]).memory / 1024 for r in results]
    random_times = [minimum(r[4]).time / 1e9 for r in results]
    random_memories = [minimum(r[4]).memory / 1024 for r in results]
    improved_errors = [r[3] for r in results]
    random_errors = [r[5] for r in results]

    p = plot(layout=(3, 1), size=(800, 1200))

    plot!(p[1], noise_levels, improved_times, label="Improved Time (s)", xlabel="Noise Level", ylabel="Time (s)", legend=:topright, title="Time Comparison (Params: $params)")
    plot!(p[1], noise_levels, random_times, label="Random Time (s)", xlabel="Noise Level", ylabel="Time (s)", legend=:topright)

    plot!(p[2], noise_levels, improved_memories, label="Improved Memory (KB)", xlabel="Noise Level", ylabel="Memory (KB)", legend=:topright, title="Memory Comparison (Params: $params)")
    plot!(p[2], noise_levels, random_memories, label="Random Memory (KB)", xlabel="Noise Level", ylabel="Memory (KB)", legend=:topright)

    plot!(p[3], noise_levels, improved_errors, label="Improved Error", xlabel="Noise Level", ylabel="Error", legend=:topright, title="Error Comparison (Params: $params)")
    plot!(p[3], noise_levels, random_errors, label="Random Error", xlabel="Noise Level", ylabel="Error", legend=:topright)

    # Save the plot to the specified path
    savefig(p, "$(output_path)_params_$(i)_noise.png")
end
