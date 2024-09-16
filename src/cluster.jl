using Clustering

"""
    struct DBSCANAlgo{T}

A structure that keeps all the necessary parameters for the DBSCAN clustering algorithm.

# Fields
- `radius::T`: The radius within which to search for neighboring points.
- `min_neighbors::Int`: The minimum number of neighbors required to form a dense region.
- `min_cluster_size::Int`: The minimum number of points required to form a cluster.

# Constructors
- `DBSCANAlgo(radius::Float64 = 0.5, min_neighbors::Int = 3, min_cluster_size::Int = 10)`: 
  Creates a new `DBSCANAlgo` instance with the specified parameters. Default values are provided for all parameters.

# Example
algo = DBSCANAlgo(radius=1.0, min_neighbors=5, min_cluster_size=15)
"""
struct DBSCANAlgo{T}
    radius::T
    min_neighbors::Int
    min_cluster_size::Int

   # Inner constructor with default values
    function DBSCANAlgo(radius::T = 0.5, min_neighbors::Int = 3, min_cluster_size::Int = 10) where T
        new{T}(radius, min_neighbors, min_cluster_size)
    end
end #struct DBSCANAlgo

function cluster_shock_points_at_t(points, dbscan_algo::DBSCANAlgo)
    dbscan_result = dbscan(points, dbscan_algo.radius, min_neighbors = dbscan_algo.min_neighbors, min_cluster_size = dbscan_algo.min_cluster_size)
    # Initial the vector to store the clusters. This vector will have the size of the number of clusters detected in all shock_positions
    shock_clusters = [] # Sequence of the cluster doesn't matter. Pushing is thread-safe
    @threads for dbscan_cluster in dbscan_result.clusters
        # Each shock_cluster is a Vector{Vector{Float64}} with n-cluster size elements
        # access to each point is done by cluster[i] which will return a vector of 2 elements [x,y]
        shock_cluster = [points[:,i] for i in dbscan_cluster.core_indices] # Only include core points to remove noises
        push!(shock_clusters, shock_cluster)
    end
    return shock_clusters
end

"""
    cluster_shock_points(dbscan_algo::DBSCANAlgo, shock_positions_over_time, flow_data)

A part of the 2D detection algorithm of function `detect`. This function clusters shock points over time using the DBSCAN algorithm.

# Arguments
- `dbscan_algo::DBSCANAlgo`: An instance of the DBSCAN clustering algorithm.
- `shock_positions_over_time`: A vector containing shock positions for each time step.
- `flow_data`: A `FlowData` object containing the flow field data.

# Returns
- `shock_clusters_over_time`: A vector containing clustered shock points for each time step.

# Description
This function clusters shock points over time using the DBSCAN algorithm. It iterates over each time step, converts shock positions to Cartesian coordinates, and applies the DBSCAN algorithm to cluster the shock points. The clustered shock points are stored in a vector, with each element corresponding to a time step.
"""
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

    @threads for t in 1:nsteps
        if isempty(shock_positions_over_time[t])
            # For frames where no shock is detected, push an empty vector
            shock_clusters_over_time[t] = []
        else
            points = cartesian_index_to_xy(shock_positions_over_time[t], x, y)
            shock_clusters_over_time[t] = cluster_shock_points_at_t(points, dbscan_algo)
        end
    end
    
    return shock_clusters_over_time
end