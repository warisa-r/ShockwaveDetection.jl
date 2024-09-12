using LsqFit
using Plots  # Import the Plots library
using BenchmarkTools  # Import the BenchmarkTools library
using Interpolations
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

    return fit, fit_error
end

# Improved function to fit a parabola with a better initial guess
function fit_parabola_with_improved_guess(xy)
    
    # Improved initial guess for parabola fitting
    x_range = LinRange(-10, 10, 100)
    y_range = LinRange(-10, 10, 100)
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

    return fit, fit_error
end

# Linear model function
function linear_model(xy, p)
    m, c = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m * x .+ c .- y
end

# Function to fit a line to the data and calculate error
function fit_line(xy)
    # Initial guess for linear fit
    p0 = [1.0, 1.0]  # Initial guesses for slope and intercept

    # Perform the fit
    fit = curve_fit(linear_model, xy, zeros(size(xy, 1)), p0)

    # Calculate error (sum of squared residuals)
    fit_error = sum(fit.resid .^ 2)

    return fit, fit_error
end

# Function to generate noisy parabolic data
function generate_noisy_parabolic_data(p)
    a, b = p
    # Generate data points along a parabola
    x = range(-10, 10, length=100)
    y = a .* x.^2 + b .* x

    # Add noise
    y_noisy = y .+ randn(length(y)) * 0.01  # Adjust noise level as needed

    # Combine into xy matrix
    xy = hcat(x, y_noisy)
    return xy
end

x_range = LinRange(-10, 10, 100)
y_range = LinRange(-10, 10, 100)
xy = generate_noisy_parabolic_data([500, 0.002])
x = xy[:, 1]
y = xy[:, 2]

#=
itp = interpolate((y_range,), y, Gridded(Linear()))
# Evaluate the gradient at the points in xy
first_derivative = only.(gradient.(Ref(itp), y_range))

b = mean(first_derivative)
println("b approximate: ", b)

# Create an interpolation object for the first derivative
itp_first_derivative = interpolate((y_range,), first_derivative, Gridded(Linear()))

# Evaluate the second derivative at the points in x
second_derivative = only.(gradient.(Ref(itp_first_derivative), y_range))
# Compute the mean of the second derivative
mean_second_derivative = mean(second_derivative)
a_estimate = mean_second_derivative/2
println("a approximate (interpolation guess): ", a_estimate)
=#

#TODO: test with LsqFit for a, b, c with the range of [0.005, 0.05, 0.1, 1 , 10, 20, 50]
#TODO: benchmark how this performs better vs random guess vs good line guess.

fit, fit_error = fit_parabola_with_improved_guess(xy)

fit_random, fit_error_random = fit_parabola_with_random_guess(xy)

println("Fit parameters (interpolation): ", fit.param)
println( "Fit error(interpolation): ", fit_error)
println("Fit parameters (random): ", fit_random.param)
println( "Fit error(random): ", fit_error_random)