module ShockwaveDetection

export write_output
export read_output_file
export convert_to_primitive
export detect_shock
export create_wave_animation, create_wave_animation_with_shock


include("dummy.jl")
include("input.jl")
include("variable_utils.jl")
include("detect_shock.jl")
include("visualize.jl")

end
