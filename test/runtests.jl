using Test
using ShockwaveDetection

@testset "Compiles" begin
    using ShockwaveDetection
end

@testset "1D shock detection" begin
    include("../examples/sod_shock_tube_1d.jl")
    rm("density_velocity_pressure_over_time_with_shock_positions.gif")
end

@testset "2D shock detection" begin
    include("../examples/sod_shock_tube_2d.jl")
end