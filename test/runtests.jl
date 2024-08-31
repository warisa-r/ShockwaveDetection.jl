using Test
using ShockwaveDetection

DATA_DIR = pkgdir(ShockwaveDetection, "examples", "data")

@testset "Compiles" begin
    using ShockwaveDetection
end

@testset "1D shock detection" begin
    flow_data = FlowData(joinpath(DATA_DIR, "sod_shock_left_1d.tape"), false)
    shock_positions_over_time = detect(flow_data, GradientShockDetectionAlgo(0.5))
    anim  = create_wave_animation_with_shock(flow_data, shock_positions_over_time)
    rm("density_velocity_pressure_over_time_with_shock_positions.gif")
end

@testset "2D shock detection" begin
    flow_data = FlowData(joinpath(DATA_DIR, "sod_shock_right_2d.tape"), false)
    point_detect_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
    dbscan_algo = DBSCANAlgo(0.25, 3, 10)
    detection = detect(flow_data, point_detect_algo, dbscan_algo)
    plot_shock_fits_over_time(flow_data, detection, true)
end