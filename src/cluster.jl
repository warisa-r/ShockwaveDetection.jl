using Clustering

# a struct that keeps all the necessary parameters for the dbscan_algo algorithm
struct DBSCANAlgo
    radius::Float64
    min_neighbors::Int
    min_cluster_size::Int

    # Inner constructor with default values
    function DBSCANAlgo(radius::Float64 = 0.5, min_neighbors::Int=3, min_cluster_size::Int=10)
        new(radius, min_neighbors, min_cluster_size)
    end
end #struct DBSCANAlgo

function cluster_shock_points_at_t(points, dbscan_algo::DBSCANAlgo)
    dbscan_result = dbscan(points, dbscan_algo.radius, min_neighbors = dbscan_algo.min_neighbors, min_cluster_size = dbscan_algo.min_cluster_size)
    # Initial the vector to store the clusters. This vector will have the size of the number of clusters detected in all shock_positions
    shock_clusters = [] # Sequence of the cluster doesn't matter. Pushing is thread-safe
    for dbscan_cluster in dbscan_result.clusters
        # Each shock_cluster is a Vector{Vector{Float64}} with n-cluster size elements
        # access to each point is done by cluster[i] which will return a vector of 2 elements [x,y]
        shock_cluster = [points[:,i] for i in dbscan_cluster.core_indices] # Only include core points to remove noises
        push!(shock_clusters, shock_cluster)
    end
    return shock_clusters
end

function cluster_shock_points(dbscan_algo::DBSCANAlgo, shock_positions_over_time, flow_data)
    bounds = flow_data.bounds
    ncells = flow_data.ncells
    nsteps = flow_data.nsteps

    # Preallocate shock_clusters_over_time with a size of nsteps, filled with undef
    shock_clusters_over_time = Vector{Any}(undef, nsteps)

    # Define the x-range
    x = range(bounds[1][1], bounds[1][2], length=ncells[1])

    # Define the y-range
    y = range(bounds[2][1], bounds[2][2], length=ncells[2])

    for t in 1:nsteps
        if isempty(shock_positions_over_time[t])
            # For framews where no shock is detected, push an empty vector
            shock_clusters_over_time[t] = []
        else
            points = cartesian_index_to_xy(shock_positions_over_time[t], x, y)
            shock_clusters_over_time[t] = cluster_shock_points_at_t(points, dbscan_algo)
        end
    end
    
    return shock_clusters_over_time
end