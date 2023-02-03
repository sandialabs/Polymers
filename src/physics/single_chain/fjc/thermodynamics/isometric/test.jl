using Test
using Polymers.Physics.SingleChain: ONE, ZERO, POINTS, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isometric: FJC

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::number_of_links" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    @test all(
        map(
            number_of_links_i ->
                FJC(
                    number_of_links_i,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                ).number_of_links,
            number_of_links_list,
        ) == number_of_links_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::link_length" begin
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * rand(parameters.number_of_loops)
    @test all(
        map(
            link_length_i ->
                FJC(
                    parameters.number_of_links_minimum,
                    link_length_i,
                    parameters.hinge_mass_reference,
                ).link_length,
            link_length_list,
        ) == link_length_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::hinge_mass" begin
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * rand(parameters.number_of_loops)
    @test all(
        map(
            hinge_mass_i ->
                FJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass_i,
                ).hinge_mass,
            hinge_mass_list,
        ) == hinge_mass_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::all_parameters" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * rand(parameters.number_of_loops)
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * rand(parameters.number_of_loops)
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                FJC(number_of_links_i, link_length_i, hinge_mass_i).number_of_links,
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) == number_of_links_list,
    )
end
