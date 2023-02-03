using Test
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isometric: FJC

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::init" begin
    @test isa(FJC(UInt8(8), 1.0, 1.0), Any)
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::number_of_links" begin
    number_of_links_list = rand(0x00:0x19, 8)
    @test all(
        map(
            number_of_links_i -> FJC(number_of_links_i, 1.0, 1.0).number_of_links,
            number_of_links_list,
        ) == number_of_links_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::link_length" begin
    link_length_list = rand(8)
    @test all(
        map(link_length_i -> FJC(0x08, link_length_i, 1.0).link_length, link_length_list) ==
        link_length_list,
    )
end
