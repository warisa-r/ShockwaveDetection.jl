using LsqFit
using Plots
using BenchmarkTools
using Interpolations

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m * x .+ b .- y
end

# Random guess line fitting function
function fit_line_with_random_guess(xy)
    # Random initial guess for line fitting
    p0 = rand(2)  # Random initial guesses for m and b
    println("Random initial guess: ", p0)

    # Perform the fit
    fit = curve_fit(line_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return fit, fit_error
end

# Improved guess line fitting function with interpolation
function fit_line_with_improved_guess(xy)
    x = xy[:, 1]
    y = xy[:, 2]

    # Interpolating for the line parameters
    itp = interpolate((x,), y, Gridded(Linear()))

    # Compute the first derivative to estimate the slope (m)
    first_derivative = only.(gradient.(Ref(itp), x))
    m_guess = mean(first_derivative)
    println("Improved initial guess for m: ", m_guess)

    # The intercept (b) is estimated as the mean of y-values
    b_guess = mean(y)
    println("Improved initial guess for b: ", b_guess)

    p0 = [m_guess, b_guess]

    # Perform the fit
    fit = curve_fit(line_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return fit, fit_error
end

# Function to generate noisy line data
function generate_noisy_line_data(m, b)
    x = range(-10, 10, length=100)
    y = m .* x .+ b
    y_noisy = y .+ randn(length(y)) * 0.01
    xy = hcat(x, y_noisy)
    return xy
end

# Separate benchmark function for line fitting with improved guess
function benchmark_line_improved_guess()
    slopes = [0.005, 0.05, 0.5, 5.0, 50.0]  # Different slopes for testing

    for m in slopes
        println("\nBenchmarking improved guess for m = $m")

        # Generate noisy data
        xy = generate_noisy_line_data(m, 0.002)  # b = 0.002 for consistency

        # Benchmark the improved guess
        @btime fit_line_with_improved_guess($xy) |> (result -> begin
            improved_fit, improved_error = result
            println("Improved Guess Parameters: ", improved_fit.param)
            println("Improved Guess Error: ", improved_error)
        end)
    end
end

# Separate benchmark function for random guess
function benchmark_line_random_guess()
    slopes = [0.005, 0.05, 0.5, 5.0, 50.0]  # Different slopes for testing

    for m in slopes
        println("\nBenchmarking random guess for m = $m")

        # Generate noisy data
        xy = generate_noisy_line_data(m, 0.002)  # b = 0.002 for consistency

        # Benchmark the random guess
        @btime fit_line_with_random_guess($xy) |> (result -> begin
            random_fit, random_error = result
            println("Random Guess Parameters: ", random_fit.param)
            println("Random Guess Error: ", random_error)
        end)
    end
end

# Run the benchmarks
benchmark_line_improved_guess()
benchmark_line_random_guess()