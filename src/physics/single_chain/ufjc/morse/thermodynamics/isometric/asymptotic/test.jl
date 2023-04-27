module Test

using Test
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ufjc.Morse.Thermodynamics.Isometric.Asymptotic: MORSEFJC

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::init" begin
    @test isa(
        MORSEFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
            parameters.link_energy_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test MORSEFJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
            parameters.link_energy_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test MORSEFJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
            parameters.link_energy_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test MORSEFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.link_stiffness_reference,
            parameters.link_energy_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::link_stiffness" begin
    for _ = 1:parameters.number_of_loops
        link_stiffness =
            parameters.link_stiffness_reference +
            parameters.link_stiffness_scale * (0.5 - rand())
        link_energy =
            parameters.link_energy_reference + parameters.link_energy_scale * (0.5 - rand())
        @test MORSEFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            link_stiffness,
            parameters.link_energy_reference,
        ).link_stiffness == link_stiffness
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::link_energy" begin
    for _ = 1:parameters.number_of_loops
        link_energy =
            parameters.link_energy_reference + parameters.link_energy_scale * (0.5 - rand())
        @test MORSEFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.link_stiffness_reference,
            link_energy,
        ).link_energy == link_energy
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isometric::asymptotic::test::base::all_parameters" begin
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
        link_energy =
            parameters.link_energy_reference + parameters.link_energy_scale * (0.5 - rand())
        link_energy =
            parameters.link_energy_reference + parameters.link_energy_scale * (0.5 - rand())
        @test all(
            MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy,
            ).number_of_links == number_of_links &&
            MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy,
            ).link_length == link_length &&
            MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy,
            ).hinge_mass == hinge_mass &&
            MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy,
            ).link_stiffness == link_stiffness &&
            MORSEFJC(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness,
                link_energy,
            ).link_energy == link_energy,
        )
    end
end

end