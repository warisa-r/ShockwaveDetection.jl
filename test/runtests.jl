using Test
using ShockwaveDetection

@testset "ShockwaveDetection.jl" begin
    @test write_output() == "Hello from Julia!"
end
