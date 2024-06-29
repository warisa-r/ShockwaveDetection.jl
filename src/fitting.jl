using LsqFit
using LinearAlgebra

struct Fitting
    model::Function
    parameters::Array{Float64}
    error::Float64
    range::Tuple{Float64, Float64}
end

# The angle of the end point and the starting point is needed to determine the proper shape
# For other models, only range of x suffices
function circle_model(xy, p)
    x0, y0, r = p
    x = xy[:, 1]
    y = xy[:, 2]
    return (x .- x0).^2 + (y .- y0).^2 .- r^2
end

# Parabola model function
function parabola_model(xy, p)
    a, b, c = p
    x = xy[:, 1]
    y = xy[:, 2]

    # Return the rotated parabola model
    return a * x.^2 .+ b * x .+ c .- y
end

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m*x .+ b .- y
end

function vline_model(xy, p)
    c = p
    x = xy[:, 1]
    return x .- c
end

function log_model(xy, p)
    a, b, c = p
    x = xy[:, 1]
    y = xy[:, 2]
    return a * log.(abs.(x .- c)) .+ b.- y
end

function fit_shock_cluster(cluster)

    function cluster_to_data_points(shock_cluster)
        len_xy = length(shock_cluster)
        xy = zeros(len_xy, 2)
        for i in 1:len_xy
            xy[i, 1] = shock_cluster[i][1] # x coordinate of the shock point in the cluster
            xy[i, 2] = shock_cluster[i][2] # y coordinate of the shock point in the cluster
        end
        return xy
    end


    xy = cluster_to_data_points(cluster)
    models = [vline_model, line_model, circle_model, parabola_model]
    #TODO: better parameter initialization from boundary conditions or information about the cluster??
    p0s = [[1.0], [1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]  # Initial parameters for each model
    
    best_fit = nothing
    least_error = Inf

    function calculate_angle(center, xy)
        x_0, y_0 = center
        x = xy[:, 1]
        y = xy[:, 2]
        return atan.(y .- y_0, x .- x_0)
    end

    for (i, model) in enumerate(models)
        p0 = p0s[i]
        fit = curve_fit(model, xy, zeros(length(cluster[:, 1])), p0)
        error = sum((fit.resid).^2)  # Calculate squared error
        
        if error < least_error
            least_error = error
            if model == circle_model
                # range of circle_model is in term of angle
                # If this model fit, assume that the cluster behave well enough for the angle estimated to be accurate
                # TODO: Is there a better way?
                x0, y0, _ = fit.param
                
                angles = calculate_angle([x0, y0], xy)
                range = (minimum(angles), maximum(angles))

                # Define a tolerance for angle closeness
                #TODO: Make this a parameter? Or fine tune this?
                tolerance = 0.1  # Adjust as needed for your application
                
                # Check if the angles span a full circle within the tolerance
                if abs(range[1]- (-pi)) < tolerance && abs(range[2] - pi) < tolerance
                    range = (pi, -pi)  # Adjust range to explicitly represent a full circle
                end
            else
                range = (minimum(xy[:, 1]), maximum(xy[:, 1]))
            end
            best_fit = Fitting(model, fit.param, error, range)
        end
    end
    print(best_fit.model)
    print(best_fit.range)

    return best_fit
end

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

function plot_shock_fit!(p, fit::Fitting, x_range)
    # TODO: Only cut from range of x_min to x_max of the cluster coordinates
    model = fit.model
    println("Model: ", fit.model)
    params = fit.parameters
    println("Parameters: ", fit.parameters)

    xy = hcat(x_range, zeros(length(x_range)))
    y_values = model(xy, params)
    plot!(p, x_range, y_values, label="Fitted line")
end

function plot_shock_fits(shock_clusters, flow_data::FlowData)
    shock_fits = fit_shock_clusters(shock_clusters)
    p = plot(title = "Shock Line Plot") # Add a title here

    bounds = flow_data.bounds
    ncells = flow_data.ncells

    x = range(bounds[1][1], stop=bounds[1][2], length=ncells[1])

    for shock_fit in shock_fits
        plot_shock_fit!(p, shock_fit, x) # Assuming x range 1:0.1:10
    end

    return p
end