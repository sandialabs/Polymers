using Polymers
using Test

@testset "Polymers.jl" begin
    @test Polymers.Physics.SingleChain.Fjc.FJC(
        UInt8(8),
        1.0,
        1.0,
    ).thermodynamics.isometric.nondimensional_force(
        0.1,
    ) > 0.0
end
