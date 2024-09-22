module NoiseAnalysis

using Distributions
using LinearAlgebra
using ShockwaveDetection #: FlowData

export NoiseData
export apply_noise, apply_noise_to_flow
export compare_shock_clusters_over_time, compare_shock_positions_over_time, compare_shock_fits_over_time

include("noise.jl")

end
