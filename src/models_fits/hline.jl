using LsqFit
using Plots
using BenchmarkTools
using Interpolations

# Horizontal line model function
function hline_model(xy, p)
    b = p[1]
    y = xy[:, 2]
    return y .- b
end

# Random guess fitting function for the horizontal line
function fit_hline_with_random_guess(xy)
    # Random initial guess for horizontal line fitting (only b)
    p0 = rand(1)  # Random initial guess for b
    println("Random initial guess: ", p0)

    # Perform the fit
    fit = curve_fit(hline_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return fit, fit_error
end

# Improved guess fitting function with interpolation for horizontal line
function fit_hline_with_improved_guess(xy)
    y = xy[:, 2]

    # Improved guess: the intercept (b) is estimated as the mean of y-values
    b_guess = mean(y)
    println("Improved initial guess for b: ", b_guess)

    p0 = [b_guess]

    # Perform the fit
    fit = curve_fit(hline_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return fit, fit_error
end

# Function to generate noisy horizontal line data
function generate_noisy_hline_data(b)
    x = range(-10, 10, length=100)
    y = fill(b, 100)  # Horizontal line with y = b
    y_noisy = y .+ randn(length(y)) * 0.01  # Adding some noise
    xy = hcat(x, y_noisy)
    return xy
end

# Separate benchmark function for horizontal line fitting with improved guess
function benchmark_hline_improved_guess()
    intercepts = [0.005, 0.05, 0.5, 5.0, 50.0]  # Different intercepts for testing

    for b in intercepts
        println("\nBenchmarking improved guess for b = $b")

        # Generate noisy data for the current intercept 'b'
        xy = generate_noisy_hline_data(b)

        # Benchmark the improved guess
        @btime fit_hline_with_improved_guess($xy) |> (result -> begin
            improved_fit, improved_error = result
            println("Improved Guess Parameter: ", improved_fit.param)
            println("Improved Guess Error: ", improved_error)
        end)
    end
end

# Separate benchmark function for random guess for horizontal line
function benchmark_hline_random_guess()
    intercepts = [0.005, 0.05, 0.5, 5.0, 50.0]  # Different intercepts for testing

    for b in intercepts
        println("\nBenchmarking random guess for b = $b")

        # Generate noisy data for the current intercept 'b'
        xy = generate_noisy_hline_data(b)

        # Benchmark the random guess
        @btime fit_hline_with_random_guess($xy) |> (result -> begin
            random_fit, random_error = result
            println("Random Guess Parameter: ", random_fit.param)
            println("Random Guess Error: ", random_error)
        end)
    end
end

# Run the benchmarks
benchmark_hline_improved_guess()
benchmark_hline_random_guess()
