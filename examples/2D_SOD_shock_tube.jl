using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/sod_shock_right_2d.tape", false)

shock_positions_over_time, angle_estimated = detect(flow_data, ImageProcessingShockDetectionAlgo(0.7, :prewitt))
create_heatmap_evo_with_shock(flow_data, shock_positions_over_time, angle_estimated, :density_field, true)