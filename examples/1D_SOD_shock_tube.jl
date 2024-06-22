using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

flow_data = FlowData("examples/data/sod_shock_orb.tape", false)

create_heatmap_evo(flow_data, :velocity_field)

