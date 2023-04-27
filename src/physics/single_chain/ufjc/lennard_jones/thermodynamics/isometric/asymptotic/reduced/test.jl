module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ufjc.LennardJones.Thermodynamics.Isometric.Asymptotic.Reduced: LENNARDJONESFJC

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::init" begin
    @test isa(
        LENNARDJONESFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test LENNARDJONESFJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test LENNARDJONESFJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test LENNARDJONESFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.link_stiffness_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::link_stiffness" begin
    for _ = 1:parameters.number_of_loops
        link_stiffness =
            parameters.link_stiffness_reference +
            parameters.link_stiffness_scale * (0.5 - rand())
        @test LENNARDJONESFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            link_stiffness,
        ).link_stiffness == link_stiffness
    end
end

@testset "physics::single_chain::ufjc::lennard_jones::thermodynamics::isometric::asymptotic::reduced::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        link_stiffness =
            parameters.link_stiffness_reference +
            parameters.link_stiffness_scale * (0.5 - rand())
        @test all(
            LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
            ).number_of_links == number_of_links &&
            LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
            ).link_length == link_length &&
            LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
            ).hinge_mass == hinge_mass &&
            LENNARDJONESFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
            ).link_stiffness == link_stiffness,
        )
    end
end

end
