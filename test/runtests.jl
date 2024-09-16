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

@testset "1D supersonic shock detection" begin
    flow_data = FlowData(joinpath(DATA_DIR, "supersonic_shock_2.tape"), false)
    point_detect_algo = GradientShockDetectionAlgo(0.2)
    detection = detect(flow_data, point_detect_algo)
    anim  = create_wave_animation_with_shock(flow_data, detection)
    rm("density_velocity_pressure_over_time_with_shock_positions.gif")
end

@testset "2D shock detection" begin
    flow_data = FlowData(joinpath(DATA_DIR, "sod_shock_right_2d.tape"), false)
    point_detect_algo = ImageProcessingShockDetectionAlgo(0.5, :prewitt)
    dbscan_algo = DBSCANAlgo(0.25, 3, 10)
    fitting_algo = FittingAlgo(0.1, false)
    detection = detect(flow_data, point_detect_algo, dbscan_algo, fitting_algo)
    plot_shock_fits_over_time(flow_data, detection, true)
end

#=
@testset "obstacle shock detection" begin
    flow_data = FlowData(joinpath(DATA_DIR, "circular_obstacle_radius_1.celltape"), false)
    point_detect_algo = ImageProcessingShockDetectionAlgo(0.2, :prewitt)
    dbscan_algo = DBSCANAlgo(0.25, 3, 10)
    detection = detect(flow_data, point_detect_algo, dbscan_algo)
    plot_shock_fits_over_time(flow_data, detection, false)
    create_heatmap_evo_with_shock(flow_data, detection, :density_field, true, false)
    rm("density_field_evolution.gif")
end
=#
