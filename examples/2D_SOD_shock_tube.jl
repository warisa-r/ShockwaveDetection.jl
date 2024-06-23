using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/sod_shock_orb.tape", false)

shock_positions_over_time, angle_estimated = detect(flow_data, ImageProcessingShockDetectionAlgo(0.5, :prewitt))
create_heatmap_evo_with_shock(flow_data, shock_positions_over_time, angle_estimated, :pressure_field)