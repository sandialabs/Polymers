module Test

using Test
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ideal.Thermodynamics.Isotensional: IDEAL

@testset "physics::single_chain::ideal::thermodynamics::isotensional::test::base::init" begin
    @test isa(
        IDEAL(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::ideal::thermodynamics::isotensional::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test IDEAL(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::ideal::thermodynamics::isotensional::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test IDEAL(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::ideal::thermodynamics::isotensional::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test IDEAL(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::ideal::thermodynamics::isotensional::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test all(
            IDEAL(number_of_links, link_length, hinge_mass).number_of_links ==
            number_of_links &&
            IDEAL(number_of_links, link_length, hinge_mass).link_length == link_length &&
            IDEAL(number_of_links, link_length, hinge_mass).hinge_mass == hinge_mass,
        )
    end
end

end
