module ShockwaveDetection

export write_output
export read_output_file, FlowData
export convert_to_primitive
export GradientRHShockDetectionAlgo
export GradientEntropyShockDetectionAlgo
export detect
export check_shock_consistency
export create_wave_animation, create_wave_animation_with_shock, create_heatmap_evo, create_heatmap_evo_with_shock


include("dummy.jl")
include("input.jl")
include("variable_utils.jl")
include("shock_detectors/1D/shock_algorithms_1D.jl") # Maybe write one files that includes all the files in the shock_detectors folder?
include("shock_analysis/normal_shock_analysis.jl")
include("visualize.jl")

end
