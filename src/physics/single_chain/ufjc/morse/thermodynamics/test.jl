module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT, PLANCK_CONSTANT
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ufjc.Morse.Thermodynamics: MORSEFJC

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::init" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::number_of_links" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::link_length" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::hinge_mass" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::link_stiffness" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::link_energy" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::base::all_parameters" begin
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

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::force" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.end_to_end_length(force, temperature)
        force_out =
            model.isometric.asymptotic.legendre.force(end_to_end_length, temperature)
        residual_abs = force - force_out
        residual_rel = residual_abs / force
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::nondimensional_force" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_force_out = model.isometric.asymptotic.legendre.nondimensional_force(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        residual_abs = nondimensional_force - nondimensional_force_out
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::helmholtz_free_energy" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.end_to_end_length(force, temperature)
        helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.legendre.helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_legendre - helmholtz_free_energy_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            )
        residual_rel = residual_abs / helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::helmholtz_free_energy_per_link" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.asymptotic.end_to_end_length_per_link(force, temperature)
        helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.gibbs_free_energy_per_link(force, temperature) +
            force * end_to_end_length_per_link
        helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.legendre.helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_per_link_legendre -
            helmholtz_free_energy_per_link_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            ) / number_of_links
        residual_rel = residual_abs / helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::relative_helmholtz_free_energy" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.end_to_end_length(force, temperature)
        relative_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.relative_gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        relative_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.legendre.relative_helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_legendre -
            relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::relative_helmholtz_free_energy_per_link" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.asymptotic.end_to_end_length_per_link(force, temperature)
        relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.relative_gibbs_free_energy_per_link(
                force,
                temperature,
            ) + force * end_to_end_length_per_link
        relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.legendre.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_per_link_legendre -
            relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::nondimensional_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length =
            model.isotensional.asymptotic.nondimensional_end_to_end_length(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.legendre.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_legendre -
            nondimensional_helmholtz_free_energy_legendre_out + (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            )
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::nondimensional_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.legendre.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link_legendre -
            nondimensional_helmholtz_free_energy_per_link_legendre_out +
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            ) / number_of_links
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length =
            model.isotensional.asymptotic.nondimensional_end_to_end_length(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_relative_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.legendre.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_legendre -
            nondimensional_relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link_legendre -
            nondimensional_relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel =
            residual_abs / nondimensional_relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::force" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.reduced.end_to_end_length(force, temperature)
        force_out = model.isometric.asymptotic.reduced.legendre.force(
            end_to_end_length,
            temperature,
        )
        residual_abs = force - force_out
        residual_rel = residual_abs / force
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::nondimensional_force" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_force_out =
            model.isometric.asymptotic.reduced.legendre.nondimensional_force(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs = nondimensional_force - nondimensional_force_out
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::helmholtz_free_energy" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.reduced.end_to_end_length(force, temperature)
        helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.reduced.gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.reduced.legendre.helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_legendre - helmholtz_free_energy_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            )
        residual_rel = residual_abs / helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::helmholtz_free_energy_per_link" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.reduced.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.end_to_end_length_per_link(
                force,
                temperature,
            )
        helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.reduced.gibbs_free_energy_per_link(
                force,
                temperature,
            ) + force * end_to_end_length_per_link
        helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.reduced.legendre.helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_per_link_legendre -
            helmholtz_free_energy_per_link_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            ) / number_of_links
        residual_rel = residual_abs / helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::relative_helmholtz_free_energy" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.reduced.end_to_end_length(force, temperature)
        relative_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.reduced.relative_gibbs_free_energy(
                force,
                temperature,
            ) + force * end_to_end_length
        relative_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.reduced.legendre.relative_helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_legendre -
            relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::relative_helmholtz_free_energy_per_link" begin
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
        end_to_end_length =
            model.isotensional.asymptotic.reduced.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.end_to_end_length_per_link(
                force,
                temperature,
            )
        relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.reduced.relative_gibbs_free_energy_per_link(
                force,
                temperature,
            ) + force * end_to_end_length_per_link
        relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.reduced.legendre.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_per_link_legendre -
            relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::nondimensional_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.reduced.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.reduced.legendre.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_legendre -
            nondimensional_helmholtz_free_energy_legendre_out + (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            )
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::nondimensional_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.reduced.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.reduced.legendre.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link_legendre -
            nondimensional_helmholtz_free_energy_per_link_legendre_out +
            (
                0.5 * log(2.0 * pi * BOLTZMANN_CONSTANT * temperature / link_stiffness) +
                log(
                    8.0 *
                    pi^2 *
                    hinge_mass *
                    link_length^2 *
                    BOLTZMANN_CONSTANT *
                    temperature / PLANCK_CONSTANT^2,
                )
            ) / number_of_links
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_legendre =
            model.isotensional.asymptotic.reduced.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_relative_helmholtz_free_energy_legendre_out =
            model.isometric.asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_legendre -
            nondimensional_relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::ufjc::morse::thermodynamics::test::legendre_asymptotic_reduced::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.asymptotic.reduced.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link_legendre -
            nondimensional_relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel =
            residual_abs / nondimensional_relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

end
