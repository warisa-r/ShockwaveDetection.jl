module ShockwaveDetection

export write_output
export read_output_file, FlowData
export convert_to_primitive, cartesian
export GradientShockDetectionAlgo
export ImageProcessingShockDetectionAlgo
export DBSCANAlgo
export cluster_shock_points
export fit_shock_clusters_over_time
export circle_model, line_model, vline_model, parabola_model, log_model
export create_wave_animation, create_wave_animation_with_shock, create_heatmap_evo, create_heatmap_evo_with_shock, plot_shock_fits_over_time
export ShockDetectionResult2D, detect


include("input.jl")
include("variable_utils.jl")
include("shock_point_detectors/1D/shock_algorithms_1D.jl") # Maybe write one files that includes all the files in the shock_detectors folder?
include("shock_point_detectors/2D/shock_algorithms_2D.jl")
include("shock_analysis/normal_shock_analysis.jl")
include("cluster.jl")
include("fitting.jl")
include("visualize.jl")
include("pipeline.jl")

end
