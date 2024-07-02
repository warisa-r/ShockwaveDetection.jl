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
    models = [vline_model, line_model, circle_model] # Use only these three firsts
    #TODO: better parameter initialization from boundary conditions or information about the cluster??
    p0s = [[1.0], [1.0, 1.0], [0.0, 0.0, 1.0]]  # Initial parameters for each model
    
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

function calculate_normal_vector(fit::Fitting, evenly_spaced_range, flow_data, t)
    density_field_t = flow_data.density_field[:,:,t]
    ncells =  flow_data.ncells
    bounds = flow_data.bounds

    if fit.model == line_model
        m, _ = fit.parameters

        # Calculate the magnitude of each vector
        magnitudes = sqrt.(1 + m ^ 2)

        # Normalize the normal vector
        # Normal vector [1, -m] or [-1, m]? 
        # TODO: Check the density behind and in front of the shock to determine the direction of the normal vector
        normals_x = 1/magnitude .* ones(length(evenly_spaced_range))
        normals_y = -m/magnitudes .* ones(length(evenly_spaced_range))
    elseif fit.model == vline_model
        # TODO: Check the density behind and in front of the shock to determine the direction of the normal vector
        c = fit.parameters[1]
        start_x = c
        # Find where c is in the x direction
        x = range(bounds[1][1], bounds[1][2], length=ncells[1])
        differences = abs.(x .- c)
        index_of_c = argmin(differences)

        # Check if the density is increasing or decreasing in the x direction
        # to identify which side is ahead of the shock (low density) and which side is behind the shock (high density)
        if density_field_t[index_of_c-5, round(Int, ncells[2] / 2)] > density_field_t[index_of_c+5, round(Int, ncells[2] / 2)]
            end_x = c + 1
        else
            end_x = c - 1
        end
        normals_x = [end_x - start_x]
        normals_y = [0]

    elseif fit.model == circle_model
        angles = evenly_spaced_range

        # Calculate normal vectors (outward from the circle center)
        normals_x = cos.(angles)
        normals_y = sin.(angles)
    end
    return normals_x, normals_y
end