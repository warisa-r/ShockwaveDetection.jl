using ShockwaveDetection
using ShockwaveProperties
using Euler2D:Euler2D

#result = Euler2D.load_cell_sim("examples/data/obstacle/circular_obstacle_radius_1.celltape")
using ImageFiltering
using Statistics: mean

# Sample matrix with NaN values
matrix_with_nan = [1.0 NaN 2.0; 3.0 4.0 5.0; 6.0 7.0 8.0]

# Function to replace NaNs with the mean of neighboring values
function replace_nan_with_mean(matrix)
    nan_replaced_matrix = copy(matrix)
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            if isnan(matrix[i, j])
                # Get the neighborhood (3x3 window) around the NaN value since the kernel size is 3x3
                imin = max(1, i-1)
                imax = min(size(matrix, 1), i+1)
                jmin = max(1, j-1)
                jmax = min(size(matrix, 2), j+1)
                
                # Extract the valid (non-NaN) neighbors
                neighbors = matrix[imin:imax, jmin:jmax]
                valid_neighbors = neighbors[.!isnan.(neighbors)]
                
                # Replace NaN with the mean of valid neighbors
                if !isempty(valid_neighbors)
                    nan_replaced_matrix[i, j] = mean(valid_neighbors)
                else # If all neighbors are NaN (in the obstacle zone)
                    nan_replaced_matrix[i, j] = 0.0
                end
            end
        end
    end
    return nan_replaced_matrix
end

# Replace NaN values in the matrix
matrix_no_nan = replace_nan_with_mean(matrix_with_nan)

# Now compute gradients using the Sobel operator
grad_dx, grad_dy = imgradients(matrix_no_nan, KernelFactors.sobel, "replicate")

# Display the gradients
display(grad_dx)
display(grad_dy)
