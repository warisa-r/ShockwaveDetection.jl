using LsqFit
using Plots  # Import the Plots library
using BenchmarkTools  # Import the BenchmarkTools library
using Interpolations
using LinearAlgebra

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m * x .+ b .- y
end

# Function to fit a line with a random initial guess
function fit_line_with_random_guess(xy)
    p0 = rand(2)  # Random initial guesses for m, b
    println("Random initial guess: ", p0)

    fit = curve_fit(line_model, xy, zeros(size(xy, 1)), p0)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Improved function to fit a line with a better initial guess
function fit_line_with_improved_guess(xy)
    x_range = LinRange(-50, 50, 100)
    y_range = LinRange(-50, 50, 100)
    x = xy[:, 1]
    y = xy[:, 2]

    c_guess = mean(y)

    itp = interpolate((x_range,), y, Gridded(Linear()))
    first_derivative = only.(gradient.(Ref(itp), y_range))

    b_guess = mean(first_derivative)
    println("b approximate: ", b_guess)

    itp_first_derivative = interpolate((y_range,), first_derivative, Gridded(Linear()))
    second_derivative = only.(gradient.(Ref(itp_first_derivative), y_range))
    mean_second_derivative = mean(second_derivative)
    a_guess = mean_second_derivative / 2
    println("a approximate (interpolation guess): ", a_guess)

    p0 = [a_guess, b_guess]

    fit = curve_fit(line_model, xy, zeros(size(xy, 1)), p0)
    fit_error = sum(fit.resid .^ 2)

    return p0, fit, fit_error
end

# Function to generate noisy linear data
function generate_noisy_linear_data(p, noise_level)
    m, b = p
    x = range(-50, 50, length=100)
    y = m .* x .+ b

    y_noisy = y .+ randn(length(y)) * noise_level
    xy = hcat(x, y_noisy)
    return xy
end

# Define the parameter combinations to test
param_combinations = [
    [0.005, 0.005],
    [1, 1],
    [500, 500],
    [0.005, 500],
    [500, 0.005],
    [0.005, 5000]
]

noise_levels = [0.0, 1.0, 5.0, 7.5, 10.0]

# Define the output path
output_path = "C:\\Users\\eivan\\Dropbox\\My PC (DESKTOP-AB2UFDR)\\Desktop\\SEP_Shockwave\\1\\"

for (i, params) in enumerate(param_combinations)
    results = []
    for noise_level in noise_levels
        xy = generate_noisy_linear_data(params, noise_level)

        improved_benchmark = @benchmark fit_line_with_improved_guess($xy)
        improved_initial_guess, improved_fit, _ = fit_line_with_improved_guess(xy)
        improved_error = norm(improved_fit.param - params)

        random_benchmark = @benchmark fit_line_with_random_guess($xy)
        random_initial_guess, random_fit, _ = fit_line_with_random_guess(xy)
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
