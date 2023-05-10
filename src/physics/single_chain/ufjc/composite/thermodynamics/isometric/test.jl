module Test

using Test
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ufjc.Composite.Thermodynamics.Isometric: CUFJC

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::init" begin
    @test isa(
        CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test CUFJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::number_of_bonds" begin
    for _ = 1:parameters.number_of_loops
        number_of_bonds =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            number_of_bonds,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).number_of_bonds == number_of_bonds
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::bond_stiffness" begin
    for _ = 1:parameters.number_of_loops
        bond_stiffness =
            parameters.bond_stiffness_reference +
            parameters.bond_stiffness_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            bond_stiffness,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).bond_stiffness == bond_stiffness
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::bond_energy" begin
    for _ = 1:parameters.number_of_loops
        bond_energy =
            parameters.bond_energy_reference + parameters.bond_energy_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            bond_energy,
            parameters.bond_scission_energy_reference,
            parameters.bond_attempt_frequency_reference,
        ).bond_energy == bond_energy
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::bond_scission_energy" begin
    for _ = 1:parameters.number_of_loops
        bond_scission_energy =
            parameters.bond_scission_energy_reference +
            parameters.bond_scission_energy_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            bond_scission_energy,
            parameters.bond_attempt_frequency_reference,
        ).bond_scission_energy == bond_scission_energy
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::bond_attempt_frequency" begin
    for _ = 1:parameters.number_of_loops
        bond_attempt_frequency =
            parameters.bond_attempt_frequency_reference +
            parameters.bond_attempt_frequency_scale * (0.5 - rand())
        @test CUFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.number_of_bonds_minimum,
            parameters.bond_stiffness_reference,
            parameters.bond_energy_reference,
            parameters.bond_scission_energy_reference,
            bond_attempt_frequency,
        ).bond_attempt_frequency == bond_attempt_frequency
    end
end

@testset "physics::single_chain::ufjc::composite::thermodynamics::isometric::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        number_of_bonds =
            rand(parameters.number_of_bonds_minimum:parameters.number_of_bonds_maximum)
        bond_stiffness =
            parameters.bond_stiffness_reference +
            parameters.bond_stiffness_scale * (0.5 - rand())
        bond_energy =
            parameters.bond_energy_reference + parameters.bond_energy_scale * (0.5 - rand())
        bond_scission_energy =
            parameters.bond_scission_energy_reference +
            parameters.bond_scission_energy_scale * (0.5 - rand())
        bond_attempt_frequency =
            parameters.bond_attempt_frequency_reference +
            parameters.bond_attempt_frequency_scale * (0.5 - rand())
        @test all(
            CUFJC(
                number_of_links,
                link_length,
                hinge_mass,
                number_of_bonds,
                bond_stiffness,
                bond_energy,
                bond_scission_energy,
                bond_attempt_frequency,
            ).number_of_links == number_of_links &&
            CUFJC(
                number_of_links,
                link_length,
                hinge_mass,
                number_of_bonds,
                bond_stiffness,
                bond_energy,
                bond_scission_energy,
                bond_attempt_frequency,
            ).link_length == link_length &&
            CUFJC(
                number_of_links,
                link_length,
                hinge_mass,
                number_of_bonds,
                bond_stiffness,
                bond_energy,
                bond_scission_energy,
                bond_attempt_frequency,
            ).hinge_mass == hinge_mass,
        )
    end
end

end
