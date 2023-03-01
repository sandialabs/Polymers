module Test

using Test
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Swfjc.Thermodynamics: SWFJC

@testset "physics::single_chain::swfjc::test::base::init" begin
    @test isa(
        SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::swfjc::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test SWFJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::swfjc::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::swfjc::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.well_width_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::swfjc::test::base::well_width" begin
    for _ = 1:parameters.number_of_loops
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            well_width,
        ).well_width == well_width
    end
end

@testset "physics::single_chain::swfjc::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        @test all(
            SWFJC(number_of_links, link_length, hinge_mass, well_width).number_of_links ==
            number_of_links &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).link_length ==
            link_length &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).hinge_mass ==
            hinge_mass &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).well_width ==
            well_width,
        )
    end
end

end
