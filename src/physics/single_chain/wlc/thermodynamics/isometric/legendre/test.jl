module Test

using Test
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Wlc.Thermodynamics.Isometric.Legendre: WLC

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::init" begin
    @test isa(
        WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test WLC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.persistance_length_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::persistance_length" begin
    for _ = 1:parameters.number_of_loops
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            persistance_length,
        ).persistance_length == persistance_length
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        @test all(
            WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length,
            ).number_of_links == number_of_links &&
            WLC(number_of_links, link_length, hinge_mass, persistance_length).link_length ==
            link_length &&
            WLC(number_of_links, link_length, hinge_mass, persistance_length).hinge_mass ==
            hinge_mass &&
            WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length,
            ).persistance_length == persistance_length,
        )
    end
end

end
