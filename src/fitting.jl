using LsqFit
using LinearAlgebra

struct Fitting
    model::Function
    parameters::Array{Float64}
    error::Float64
    range::Tuple{Float64, Float64}
end

#TODO: hline_modelsh

# Circle model function
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
    return a * x.^2 .+ b * x .+ c .- y
end

# Line model function
function line_model(xy, p)
    m, b = p
    x = xy[:, 1]
    y = xy[:, 2]
    return m*x .+ b .- y
end

# Horizontal line model function
function hline_model(xy, p)
    b = p
    y = xy[:, 2]
    return y .- b
end

# Vertical line model function
function vline_model(xy, p)
    c = p
    x = xy[:, 1]
    return x .- c
end

# Helper function to initialize better initial guess based on boundary conditions or prior knowledge
function compute_better_initial_guess(cluster)
    # Add logic to compute a better initial guess, for now we use predefined values
    p0s = [[1.0], [1.0], [1.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
    return p0s
end

# Helper function to generate random initial guesses
function generate_random_initial_guess(cluster)
    p0s = [[rand()], [rand()], [rand(), rand()], [rand(), rand(), rand()], [rand(), rand(), rand()]]
    return p0s
end

function fit_shock_cluster(cluster, make_initial_guess::Bool)
    # Helper function to convert a cluster to a matrix of data points
    function cluster_to_data_points(shock_cluster)
        len_xy = length(shock_cluster)
        xy = zeros(len_xy, 2)
        @threads for i in 1:len_xy
            xy[i, 1] = shock_cluster[i][1]
            xy[i, 2] = shock_cluster[i][2]
        end
        return xy
    end

    xy = cluster_to_data_points(cluster)
    models = [vline_model, hline_model, line_model, circle_model, parabola_model]
    
    # Select initial parameters based on whether we use better initial guess or random guess
    p0s = make_initial_guess ? compute_better_initial_guess(cluster) : generate_random_initial_guess(cluster)
    
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
        error = sum((fit.resid).^2)
        
        if error < least_error
            least_error = error
            if model == circle_model
                x0, y0, _ = fit.param
                angles = calculate_angle([x0, y0], xy)
                range = (minimum(angles), maximum(angles))
                tolerance = 0.1
                
                if abs(range[1]- (-pi)) < tolerance && abs(range[2] - pi) < tolerance
                    range = (pi, -pi)
                end
            else
                range = (minimum(xy[:, 1]), maximum(xy[:, 1]))
            end
            best_fit = Fitting(model, fit.param, error, range)
        end
    end

    return best_fit
end

function fit_shock_clusters(shock_clusters::Vector{Any}, make_initial_guess::Bool)
    shock_fits = []
    if !isempty(shock_clusters)
        @threads for shock_cluster in shock_clusters
            best_fit = fit_shock_cluster(shock_cluster, make_initial_guess)
            push!(shock_fits, best_fit)
        end
    end
    return shock_fits
end

function fit_shock_clusters_over_time(shock_clusters_over_time::Vector{Any}, make_initial_guess::Bool)
    nsteps = length(shock_clusters_over_time)
    shock_fits_over_time = Vector{Any}(undef, nsteps)
    @threads for t in 1:nsteps
        shock_fits = []
        if !isempty(shock_clusters_over_time[t])
            shock_fits = fit_shock_clusters(shock_clusters_over_time[t], make_initial_guess)
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
        m, b = fit.parameters
        m = m == 0 ? 1e-6 : m
        magnitude = sqrt.(1 + m^2)
        normal_quantity_x = 1/magnitude
        normal_quantity_y = m/magnitude
        x = range(bounds[1][1], bounds[1][2], length=ncells[1])
        y = range(bounds[2][1], bounds[2][2], length=ncells[2])
        mid_x = evenly_spaced_range[round(Int, length(evenly_spaced_range) / 2)]
        mid_y = m * mid_x + b
        x_differences = abs.(x .- mid_x)
        y_differences = abs.(y .- mid_y)
        index_of_mid_x = argmin(x_differences)
        index_of_mid_y = argmin(y_differences)
        coordinate_movement_x = m < 1 ? round(Int, 1/m) : 1
        coordinate_movement_y = m < 1 ? 1 : round(Int, m)
        density_left = density_field_t[index_of_mid_x - coordinate_movement_x, index_of_mid_y + coordinate_movement_y]
        density_right = density_field_t[index_of_mid_x + coordinate_movement_x, index_of_mid_y - coordinate_movement_y]
        if density_left > density_right
            normal_quantity_y = -normal_quantity_y
        else
            normal_quantity_x = -normal_quantity_x
        end
        normals_x = normal_quantity_x .* ones(length(evenly_spaced_range))
        normals_y = normal_quantity_y .* ones(length(evenly_spaced_range))
    elseif fit.model == vline_model
        c = fit.parameters[1]
        start_x = c
        x = range(bounds[1][1], bounds[1][2], length=ncells[1])
        differences = abs.(x .- c)
        index_of_c = argmin(differences)
        end_x = density_field_t[index_of_c-5, round(Int, ncells[2] / 2)] > density_field_t[index_of_c+5, round(Int, ncells[2] / 2)] ? c + 1 : c - 1
        normals_x = [end_x - start_x]
        normals_y = [0]
    elseif fit.model == hline_model
        b = fit.parameters[1]
        start_y = b
        y = range(bounds[2][1], bounds[2][2], length=ncells[2])
        differences = abs.(y .- b)
        index_of_b = argmin(differences)
        end_y = density_field_t[round(Int, ncells[1] / 2), index_of_b-5] > density_field_t[round(Int, ncells[1] / 2), index_of_b+5] ? b + 1 : b - 1
        normals_x = [0]
        normals_y = [end_y - start_y]
    elseif fit.model == circle_model
        angles = evenly_spaced_range
        normals_x = cos.(angles)
        normals_y = sin.(angles)
    elseif fit.model == parabola_model
        a, b, _ = fit.parameters
        function derivative_parabola(x)
            return 2 * a * x + b
        end
        normals_x = similar(evenly_spaced_range)
        normals_y = similar(evenly_spaced_range)
        for (i, x) in enumerate(evenly_spaced_range)
            slope = derivative_parabola(x)
            if slope != 0
                slope_normal = -1/slope
                if a > 0
                    normals_x[i] = 1/sqrt(1 + slope_normal^2)
                    normals_y[i] = slope_normal/sqrt(1 + slope_normal^2)
                    if x < 0
                        normals_x[i] = -normals_x[i]
                        normals_y[i] = -normals_y[i]
                    end
                else
                    normals_x[i] = -1/sqrt(1 + slope_normal^2)
                    normals_y[i] = -slope_normal/sqrt(1 + slope_normal^2)
                    if x > 0
                        normals_x[i] = -normals_x[i]
                        normals_y[i] = -normals_y[i]
                    end
                end
            else
                normals_x[i] = 0
                normals_y[i] = a > 0 ? -1 : 1
            end
        end
    end
    return normals_x, normals_y
end
