module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ONE, ZERO, POINTS, integrate, parameters
using Polymers.Physics.SingleChain.Ufjc.Morse.Thermodynamics.Isotensional.Legendre: MORSEFJC

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::init" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::number_of_links" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::link_length" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::hinge_mass" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::link_stiffness" begin
    for _ = 1:parameters.number_of_loops
        link_stiffness =
            parameters.link_stiffness_reference +
            parameters.link_stiffness_scale * (0.5 - rand())
        @test MORSEFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            link_stiffness,
            parameters.link_energy_reference,
        ).link_stiffness == link_stiffness
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::link_energy" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::base::all_parameters" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::nondimensional::helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_helmholtz_free_energy =
            model.nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        helmholtz_free_energy = model.helmholtz_free_energy(force, temperature)
        residual_abs =
            helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::nondimensional::helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(force, temperature)
        residual_abs =
            helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::nondimensional::relative_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_force,
                temperature,
            )
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(force, temperature)
        residual_abs =
            relative_helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::nondimensional::relative_helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(force, temperature)
        residual_abs =
            relative_helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::per_link::helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        helmholtz_free_energy = model.helmholtz_free_energy(force, temperature)
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(force, temperature)
        residual_abs =
            helmholtz_free_energy / number_of_links - helmholtz_free_energy_per_link
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::per_link::relative_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(force, temperature)
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(force, temperature)
        residual_abs =
            relative_helmholtz_free_energy / number_of_links -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::per_link::nondimensional_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_helmholtz_free_energy =
            model.nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_force,
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::per_link::nondimensional_relative_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_force,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy / number_of_links -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::relative::helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        helmholtz_free_energy = model.helmholtz_free_energy(force, temperature)
        helmholtz_free_energy_0 = model.helmholtz_free_energy(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(force, temperature)
        residual_abs =
            helmholtz_free_energy - helmholtz_free_energy_0 - relative_helmholtz_free_energy
        residual_rel = residual_abs / helmholtz_free_energy_0
        @test abs(residual_abs) <=
              BOLTZMANN_CONSTANT * temperature * number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::relative::helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(force, temperature)
        helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(force, temperature)
        residual_abs =
            helmholtz_free_energy_per_link - helmholtz_free_energy_per_link_0 -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / helmholtz_free_energy_per_link_0
        @test abs(residual_abs) <= BOLTZMANN_CONSTANT * temperature * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::relative::nondimensional_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_helmholtz_free_energy =
            model.nondimensional_helmholtz_free_energy(nondimensional_force, temperature)
        nondimensional_helmholtz_free_energy_0 =
            model.nondimensional_helmholtz_free_energy(ZERO, temperature)
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_force,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy - nondimensional_helmholtz_free_energy_0 -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_0
        @test abs(residual_abs) <= number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::relative::nondimensional_helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force_max =
            sqrt(link_stiffness * link_energy / 8) / BOLTZMANN_CONSTANT / temperature *
            link_length
        nondimensional_force = nondimensional_force_max * rand()
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_0 =
            model.nondimensional_helmholtz_free_energy_per_link(ZERO, temperature)
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link -
            nondimensional_helmholtz_free_energy_per_link_0 -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_0
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::zero::relative_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        @test abs(relative_helmholtz_free_energy_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature * number_of_links
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::zero::relative_helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_helmholtz_free_energy_per_link_0 =
            model.relative_helmholtz_free_energy_per_link(
                ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
                temperature,
            )
        @test abs(relative_helmholtz_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::zero::nondimensional_relative_helmholtz_free_energy" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_0 =
            model.nondimensional_relative_helmholtz_free_energy(ZERO, temperature)
        @test abs(nondimensional_relative_helmholtz_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::isotensional::legendre::test::zero::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
        model =
            MORSEFJC(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_per_link_0 =
            model.nondimensional_relative_helmholtz_free_energy_per_link(ZERO, temperature)
        @test abs(nondimensional_relative_helmholtz_free_energy_per_link_0) <= ZERO
    end
end

end
