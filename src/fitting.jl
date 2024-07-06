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
function initialize_circle_params(xy)
    x0 = mean(xy[:, 1])
    y0 = mean(xy[:, 2])
    r = mean(sqrt.((xy[:, 1] .- x0).^2 .+ (xy[:, 2] .- y0).^2))
    return [x0, y0, r]
end

# Parabola model function
function initialize_parabola_params(xy)
    x = xy[:, 1]
    y = xy[:, 2]
    # Return the rotated parabola model
    p = polyfit(x,y,2) #Fit a second-order polynomial    
    return p.coeffs
end

# Line model function
function initialize_line_params(xy)
    x = xy[:, 1]
    y = xy[:, 2]
    A = [x ones(length(x))]
    β = A \ y
    return β
end

function initialize_vline_params(xy)
    c = mean(xy[:,1])
    return [c]
end

function initialize_log_params(xy)
    x = xy[:, 1]
    y = xy[:, 2]
    a = (maximum(y)- minimum(y)) / (log(maximum(x)) - log(minimum(x) + 1e-6))
    b = mean(y) 
    c = mean(x)
    return [a,b,c]
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
    p0s = [initialize_vline_params(xy), initialize_line_params(xy),
    initialize_circle_params(xy), initialize_parabola_params(xy),
    initialize_log_params(xy)
    ]  
    # Initial parameters for each model
    
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