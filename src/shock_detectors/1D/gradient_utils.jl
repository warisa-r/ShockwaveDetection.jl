using Interpolations: interpolate, Gridded, Linear, gradient

# Function to compute gradients using Interpolations.jl
# without explicit function definition
function compute_gradients(arr::AbstractVector, x)
    # Convert StepRangeLen to LinRange because that's what Interpolations.jl expects
    x_linrange = LinRange(first(x), last(x), length(x))
    itp = interpolate((x_linrange,), arr, Gridded(Linear()))
    grad = only.(gradient.(Ref(itp), x_linrange))
    return grad
end