struct EntropyGradientShockDetectionAlgo <: Abstract1DShockDetectionAlgo
    threshold::Float64
    entropy::Float64
end # Simple1DShockDetector