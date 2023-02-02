using Test
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isometric: FJC

@testset "physics::single_chain::thermodynamics::isometric::base::init" begin
    for _ = 0:7
        @test isa(FJC(UInt8(8), 1.0, 1.0), Any)
    end
end
@testset "physics::single_chain::thermodynamics::isometric::base::number_of_links" begin
    for number_of_links in rand(0x00:0x19, 8)
        @test FJC(number_of_links, 1.0, 1.0).number_of_links == number_of_links
    end
end
