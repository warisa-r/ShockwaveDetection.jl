using LsqFit
using LinearAlgebra
using Plots

struct Fitting
    model::Function
    parameters::Array{Float64}
    error::Float64
end

function circle_model(xy, p)
    x0, y0, r = p
    x = xy[:, 1]
    y = xy[:, 2]
    return (x .- x0).^2 + (y .- y0).^2 .- r^2
end

# Parabola model function
function parabola_model(xy, p)
    a, x_0, y_0, theta = p
    x = xy[:, 1]
    y = xy[:, 2]

    # Return the rotated parabola model
    return -x.*sin(theta)+ y.*cos(theta) - a .* (x.*cos(theta)+y.*sin(theta) .- x_0).^2 .+ y_0
end

function ellipse_model(xy, p)
    x0, y0, a, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return (x .- x0).^2 ./a.^2 + (y .- y0).^2 ./b.^2 .- 1
end

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m*x .+ b .- y
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
    models = [line_model, circle_model]
    #TODO: better parameter initialization from boundary conditions or information about the cluster??
    p0s = [[1.0, 1.0], [0.0, 0.0, 1.0]]  # Initial parameters for each model
    
    best_fit = nothing
    least_error = Inf

    for (i, model) in enumerate(models)
        p0 = p0s[i]
        fit = curve_fit(model, xy, zeros(length(cluster[:, 1])), p0)
        error = sum((fit.resid).^2)  # Calculate squared error
        
        if error < least_error
            least_error = error
            best_fit = Fitting(model, fit.param, error)
        end
    end

    return best_fit
end

function fit_shock_clusters(shock_clusters)
    shock_lines = []
    if !isempty(shock_clusters)
        for shock_cluster in shock_clusters
            best_fit = fit_shock_cluster(shock_cluster)
            push!(shock_lines, best_fit)
        end
    end
    return shock_lines
end

function fit_shock_clusters_over_time(shock_clusters_over_time)
    shock_lines_over_time = []
    for shock_clusters in shock_clusters_over_time
        shock_lines = fit_shock_clusters(shock_clusters)
        push!(shock_lines_over_time, shock_lines)
    end
    return shock_lines_over_time 
end

function plot_shock_line!(p, fit::Fitting, x_range)
    model = fit.model
    println("Model: ", fit.model)
    params = fit.parameters
    println("Parameters: ", fit.parameters)

    xy = hcat(x_range, zeros(length(x_range)))
    y_values = model(xy, params)
    plot!(p, x_range, y_values, label="Fitted line")
end

function plot_shock_lines(shock_clusters, flow_data::FlowData)
    shock_lines = fit_shock_clusters(shock_clusters)
    p = plot(title = "Shock Line Plot") # Add a title here

    bounds = flow_data.bounds
    ncells = flow_data.ncells

    x = range(bounds[1][1], stop=bounds[1][2], length=ncells[1])

    for shock_line in shock_lines
        plot_shock_line!(p, shock_line, x) # Assuming x range 1:0.1:10
    end

    display(p)
    
end