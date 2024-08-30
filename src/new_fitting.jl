using Statistics
using LinearAlgebra
using LsqFit
using Optim

struct Fitting
    model::Function
    parameters::Array{Float64}
    error::Float64
    range::Tuple{Float64, Float64}
end

# Function to compute the geometric median
function geometric_median(xy, tol=1e-5, max_iter=100)
    median = mean(xy, dims=1)
    for _ in 1:max_iter
        distances = sqrt.(sum((xy .- median).^2, dims=2))
        non_zero = distances .> 0
        if sum(non_zero) == 0
            return median
        end
        inv_distances = 1.0 ./ distances[non_zero]
        weighted_sum = sum(xy[non_zero, :] .* inv_distances, dims=1) / sum(inv_distances)
        shift = weighted_sum - median
        median += shift
        if norm(shift) < tol
            break
        end
    end
    return median
end

# Function to initialize parameters for the horizontal line model
function initialize_hline_params(xy)
    c = mean(xy[:, 2])  # Initialize c to the mean of the y-coordinates
    return [c]
end

# Improved parameter initialization function for the vertical line model
function initialize_vline_parameters(xy)
    x = xy[:, 1]
    
    # Estimate c as the mean or median of x values
    c_initial = mean(x)  # Or use median(x) if the data is skewed
    
    return [c_initial]  # Return as an array for consistency with other models
end

# Vertical line model function
function vline_model(xy, p)
    c = p[1]  # Parameter c is the x-coordinate of the vertical line
    x = xy[:, 1]  # Extract the x-coordinates of the data points
    return x .- c  # Residuals between the model and the actual data
end

# Objective function to minimize for vertical line fitting
function objective(p, xy)
    residuals = vline_model(xy, p)
    return sum(residuals.^2)  # Sum of squared residuals
end

# Updated function to initialize parameters for the line model
function initialize_line_params(xy)
    x = xy[:, 1]
    y = xy[:, 2]

    # Create the design matrix A, where the first column is x, and the second column is ones
    A = [x ones(length(x))]

    # Solve for the parameters β using the backslash operator
    β = A \ y

    # β contains [m, b] where m is the slope and b is the intercept
    m = β[1]
    b = β[2]

    return [m, b]
end

# Function to initialize parameters for the circle model
function initialize_circle_params(xy)
    center = geometric_median(xy)
    x0, y0 = center[1], center[2]
    distances = sqrt.((xy[:,1] .- x0).^2 .+ (xy[:,2] .- y0).^2)
    r_initial = median(distances)
    return [x0, y0, r_initial]
end

# Function to initialize parameters for the parabola model
function initialize_parabola_params(xy)
    x = xy[:, 1]
    y = xy[:, 2]

    # Initial guesses based on vertex form
    h = mean(x)
    k = mean(y)

    # Estimate 'a' by considering the spread of y-values around the vertex
    # Use the difference between max and min y-values relative to mean x
    a_initial = (maximum(y) - minimum(y)) / (maximum(x) - minimum(x))^2

    # Convert vertex form to standard form: y = ax^2 + bx + c
    a = a_initial
    b = -2 * a * h
    c = a * h^2 + k

    return [a, b, c]
end

# Function to convert a cluster to a matrix of data points
function cluster_to_data_points(shock_cluster)
    len_xy = length(shock_cluster)
    xy = zeros(len_xy, 2)
    for i in 1:len_xy
        xy[i, 1] = shock_cluster[i][1]  # x coordinate
        xy[i, 2] = shock_cluster[i][2]  # y coordinate
    end
    return xy
end

# Horizontal line model function
function hline_model(xy, p)
    c = p[1]  # y-coordinate of the horizontal line
    y = xy[:, 2]  # y-coordinates of the data points
    return c .- y  # Residuals between the model and the actual data
end

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m * x .+ b .- y  # Residuals between the model and the actual data
end

# Circle model function
function circle_model(xy, p)
    x0, y0, r = p
    x = xy[:, 1]
    y = xy[:, 2]
    return (x .- x0).^2 + (y .- y0).^2 .- r^2  # Residuals for the circle model
end

# Parabola model function
function parabola_model(xy, p)
    a, b, c = p
    x = xy[:, 1]
    y = xy[:, 2]
    return a * x.^2 .+ b * x .+ c .- y  # Residuals for the parabola model
end

# Function to calculate angles, useful for circular fits
function calculate_angle(center, xy)
    x_0, y_0 = center
    x = xy[:, 1]
    y = xy[:, 2]
    return atan.(y .- y_0, x .- x_0)
end

# Function to fit a single shock cluster
function fit_shock_cluster(cluster)
    xy = cluster_to_data_points(cluster)
    models = [hline_model, vline_model, line_model, circle_model, parabola_model]

    # Initial parameter guesses for each model
    p0s = [
        initialize_hline_params(xy),
        initialize_vline_parameters(xy),
        initialize_line_params(xy),
        initialize_circle_params(xy),
        initialize_parabola_params(xy)
    ]

    best_fit = nothing
    least_error = Inf

    for (i, model) in enumerate(models)
        p0 = p0s[i]
        fit = curve_fit(model, xy, zeros(length(cluster[:, 1])), p0)
        error = sum((fit.resid).^2)  # Calculate squared error

        if error < least_error
            least_error = error
            if model == circle_model
                x0, y0, _ = fit.param
                angles = calculate_angle([x0, y0], xy)
                range = (minimum(angles), maximum(angles))

                tolerance = 0.1  # Adjust as needed

                if abs(range[1] - (-π)) < tolerance && abs(range[2] - π) < tolerance
                    range = (π, -π)  # Adjust range for a full circle
                end
            else
                range = (minimum(xy[:, 1]), maximum(xy[:, 1]))
            end

            best_fit = Fitting(model, fit.param, error, range)
        end
    end

    return best_fit
end

# Function to fit shock clusters
function fit_shock_clusters(shock_clusters)
    shock_fits = []
    if !isempty(shock_clusters)
        for shock_cluster in shock_clusters
            best_fit = fit_shock_cluster(shock_cluster)
            push!(shock_fits, best_fit)
        end
    end
    return shock_fits
end

# Function to fit shock clusters over time
function fit_shock_clusters_over_time(shock_clusters_over_time)
    shock_fits_over_time = []
    for shock_clusters in shock_clusters_over_time
        shock_fits = []
        if !isempty(shock_clusters)
            shock_fits = fit_shock_clusters(shock_clusters)
        end
        push!(shock_fits_over_time, shock_fits)
    end
    return shock_fits_over_time
end
