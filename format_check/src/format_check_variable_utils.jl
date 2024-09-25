"""
    cartesian_index_to_xy(shock_positions_t, x, y) -> Bool

Converts Cartesian indices of shock positions into their corresponding x and y coordinates.

# Arguments
- `shock_positions_t`: An array of CartesianIndex indicating the positions of shocks.
- `x`: An array representing the x-coordinates in the domain.
- `y`: An array representing the y-coordinates in the domain.

# Returns
- `Bool`: Returns `true` if the conversion is successful, otherwise `false`.
"""
function cartesian_index_to_xy(shock_positions_t, x, y)
    if isempty(shock_positions_t) || isempty(x) || isempty(y)
        return false  # Return false if input arrays are empty
    end

    coordinate_matrix = zeros(eltype(x), 2, length(shock_positions_t))
    
    # Extract x and y coordinates from CartesianIndex
    xs = [pos[1] for pos in shock_positions_t]
    ys = [pos[2] for pos in shock_positions_t]

    # Check if indices are valid
    if any(x_pos < 1 || x_pos > length(x) || y_pos < 1 || y_pos > length(y) for (x_pos, y_pos) in zip(xs, ys))
        return false  # Return false if any index is out of bounds
    end

    # Nail down where this is in x and y before scattering since shock_positions is just sets of indices not actual x and y values
    x_shocks = [x[x_pos] for x_pos in xs]
    y_shocks = [y[y_pos] for y_pos in ys]

    # Assign x_shocks to the first row and y_shocks to the second row of coordinate_matrix
    coordinate_matrix[1, :] = x_shocks
    coordinate_matrix[2, :] = y_shocks

    return true  # Return true if the conversion was successful
end

# Sample data
shock_positions_t = CartesianIndex[CartesianIndex(1, 1), CartesianIndex(2, 2)]
x = [0.0, 1.0, 2.0, 3.0]
y = [0.0, 1.0, 2.0, 3.0]

# Call the function
result = cartesian_index_to_xy(shock_positions_t, x, y)

# Print the result
println("Conversion successful: ", result)  # This will print true or false
