module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ZERO, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.ModifiedCanonical.Asymptotic.StrongPotential:
    FJC

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test FJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test FJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test all(
            FJC(number_of_links, link_length, hinge_mass).number_of_links ==
            number_of_links &&
            FJC(number_of_links, link_length, hinge_mass).link_length == link_length &&
            FJC(number_of_links, link_length, hinge_mass).hinge_mass == hinge_mass,
        )
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::nondimensional::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force = model.nondimensional_force(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
        )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        force = model.force(potential_distance, potential_stiffness, temperature)
        residual_abs =
            force / BOLTZMANN_CONSTANT / temperature * link_length - nondimensional_force
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::nondimensional::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        helmholtz_free_energy = model.helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::nondimensional::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::nondimensional::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            relative_helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::nondimensional::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(
                potential_distance,
                potential_stiffness,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::per_link::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        helmholtz_free_energy = model.helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            helmholtz_free_energy / number_of_links - helmholtz_free_energy_per_link
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::per_link::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(
                potential_distance,
                potential_stiffness,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy / number_of_links -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::per_link::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy / number_of_links -
            nondimensional_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::per_link::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy / number_of_links -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::relative::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        helmholtz_free_energy = model.helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        helmholtz_free_energy_0 = model.helmholtz_free_energy(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            helmholtz_free_energy - helmholtz_free_energy_0 - relative_helmholtz_free_energy
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::relative::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(
                potential_distance,
                potential_stiffness,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_per_link - helmholtz_free_energy_per_link_0 -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::relative::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(
            ZERO,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy - nondimensional_helmholtz_free_energy_0 -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::relative::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small * (1.0 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_0 =
            model.nondimensional_helmholtz_free_energy_per_link(
                ZERO,
                nondimensional_potential_stiffness,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link -
            nondimensional_helmholtz_free_energy_per_link_0 -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::zero::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        @test abs(relative_helmholtz_free_energy_0) <=
              ZERO * number_of_links * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::zero::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_helmholtz_free_energy_per_link_0 =
            model.relative_helmholtz_free_energy_per_link(
                ZERO * number_of_links * link_length,
                potential_stiffness,
                temperature,
            )
        @test abs(relative_helmholtz_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::zero::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_0 =
            model.nondimensional_relative_helmholtz_free_energy(
                ZERO,
                nondimensional_potential_stiffness,
            )
        @test abs(nondimensional_relative_helmholtz_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::zero::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_per_link_0 =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                ZERO,
                nondimensional_potential_stiffness,
            )
        @test abs(nondimensional_relative_helmholtz_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::connection::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / (number_of_links * link_length)^2 *
            BOLTZMANN_CONSTANT *
            temperature
        force = model.force(potential_distance, potential_stiffness, temperature)
        h = parameters.rel_tol * number_of_links * link_length
        force_from_derivative =
            (
                model.relative_helmholtz_free_energy(
                    potential_distance + 0.5 * h,
                    potential_stiffness,
                    temperature,
                ) - model.relative_helmholtz_free_energy(
                    potential_distance - 0.5 * h,
                    potential_stiffness,
                    temperature,
                )
            ) / h
        residual_abs = force - force_from_derivative
        residual_rel = residual_abs / force
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::test::connection::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_small
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force = model.nondimensional_force(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
        )
        h = parameters.rel_tol
        nondimensional_force_from_derivative =
            (
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_potential_distance + 0.5 * h,
                    nondimensional_potential_stiffness,
                ) - model.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_potential_distance - 0.5 * h,
                    nondimensional_potential_stiffness,
                )
            ) / h
        residual_abs = nondimensional_force - nondimensional_force_from_derivative
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= h
    end
end

end
