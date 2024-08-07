using LsqFit
using LinearAlgebra

struct Fitting
    model::Function
    parameters::Array{Float64}
    error::Float64
    range::Tuple{Float64, Float64}
end

#TODO: hline_model

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

function fit_shock_cluster(cluster)

    # Helper function to convert a cluster to a matrix of data points in a form that LsqFit can use
    function cluster_to_data_points(shock_cluster)
        len_xy = length(shock_cluster)
        xy = zeros(len_xy, 2)
        @threads for i in 1:len_xy
            xy[i, 1] = shock_cluster[i][1] # x coordinate of the shock point in the cluster
            xy[i, 2] = shock_cluster[i][2] # y coordinate of the shock point in the cluster
        end
        return xy
    end


    xy = cluster_to_data_points(cluster)
    models = [vline_model, line_model, circle_model, parabola_model] # Use only these three firsts
    #TODO: better parameter initialization from boundary conditions or information about the cluster??
    p0s = [[1.0], [1.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]  # Initial parameters for each model
    
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
        @threads for shock_cluster in shock_clusters # Sequence of fit doesn't matter
            best_fit = fit_shock_cluster(shock_cluster)
            push!(shock_fits, best_fit)
        end
    end
    return shock_fits
end

function fit_shock_clusters_over_time(shock_clusters_over_time)
    nsteps = length(shock_clusters_over_time)
    shock_fits_over_time = Vector{Any}(undef, nsteps)
    @threads for t in 1:nsteps
        shock_fits = []
        if !isempty(shock_clusters_over_time[t])
            shock_fits = fit_shock_clusters(shock_clusters_over_time[t])
        end
        shock_fits_over_time[t] = shock_fits
    end
    return shock_fits_over_time 
end

function calculate_normal_vector(fit::Fitting, evenly_spaced_range, flow_data, t)
    density_field_t = flow_data.density_field[:,:,t]
    ncells =  flow_data.ncells
    bounds = flow_data.bounds

    if fit.model == line_model
        #TODO: Check if this work! This is just pure implementation no test!
        m, b = fit.parameters

        if m == 0
            m = 1e-6
            #TODO: hline model
        end

        # Calculate the magnitude of each vector
        magnitude = sqrt.(1 + m ^ 2)

        normal_quantity_x = 1/magnitude
        normal_quantity_y = m/magnitude

        x = range(bounds[1][1], bounds[1][2], length=ncells[1])
        y = range(bounds[2][1], bounds[2][2], length=ncells[2])
        
        # Find the point in the middle of the line
        mid_x = evenly_spaced_range[round(Int, length(evenly_spaced_range) / 2)]
        mid_y = m * mid_x + b

        # Find out the index of these points in the x and y direction
        x_differences = abs.(x .- mid_x)
        y_differences = abs.(y .- mid_y)

        index_of_mid_x = argmin(x_differences)
        index_of_mid_y = argmin(y_differences)

        if m < 1
            coordinate_movement_x = round(Int, 1/m)
            coordinate_movement_y = 1
        else
            coordinate_movement_x = 1
            coordinate_movement_y = round(Int, m)
        end

        density_left = density_field_t[index_of_mid_x - coordinate_movement_x, index_of_mid_y + coordinate_movement_y]
        density_right = density_field_t[index_of_mid_x + coordinate_movement_x, index_of_mid_y - coordinate_movement_y]

        # Shock is moving from left to right
        if density_left > density_right
            normal_quantity_x = normal_quantity_x
            normal_quantity_y = -normal_quantity_y
        else
            normal_quantity_x = -normal_quantity_x
            normal_quantity_y = normal_quantity_y
        end

        normals_x = normal_quantity_x .* ones(length(evenly_spaced_range))
        normals_y = normal_quantity_y .* ones(length(evenly_spaced_range))
    elseif fit.model == vline_model
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
    elseif fit.model == parabola_model
        a, b, _ = fit.parameters
        function derivative_parabola(x)
            return 2 * a * x + b
        end

        # Initialize the array of normal vectors
        normals_x = similar(evenly_spaced_range)
        normals_y = similar(evenly_spaced_range)
        
        for (i, x) in enumerate(evenly_spaced_range)
            # Calculate the slope of the tangent line at the point
            slope = derivative_parabola(x)

            if slope != 0
                slope_normal = -1/slope
                if a > 0
                    normals_x[i] = 1/sqrt(1 + slope_normal^2)
                    normals_y[i] = slope_normal/sqrt(1 + slope_normal^2)
                    if x < 0
                        # Adjust direction if on the left side of the vertex
                        normals_x[i] = -normals_x[i]
                        normals_y[i] = -normals_y[i]
                    end
                else
                    normals_x[i] = -1/sqrt(1 + slope_normal^2)
                    normals_y[i] = -slope_normal/sqrt(1 + slope_normal^2)
                    if x > 0
                        # Adjust direction if on the right side of the vertex
                        normals_x[i] = -normals_x[i]
                        normals_y[i] = -normals_y[i]
                    end
                end
            else
                # If the slope is zero, the normal vector is vertical
                if a > 0
                    normals_x[i] = 0
                    normals_y[i] = -1
                else
                    normals_x[i] = 0
                    normals_y[i] = 1
                end
            end
        end
    end
    return normals_x, normals_y
end