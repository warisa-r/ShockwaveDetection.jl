
# Abstract base type for 1D shock detectors
abstract type Abstract1DShockDetectionAlgo end

include("gradient_utils.jl")
include("gradient_rhs_1D.jl")
include("gradient_entropy_1D.jl")