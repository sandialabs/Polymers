module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ZERO, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.ModifiedCanonical.Asymptotic.WeakPotential:
    FJC

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::base::number_of_links" begin
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

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::base::link_length" begin
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

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::base::hinge_mass" begin
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

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::base::all_parameters" begin
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

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
        )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        end_to_end_length =
            model.end_to_end_length(potential_distance, potential_stiffness, temperature)
        residual_abs = end_to_end_length / link_length - nondimensional_end_to_end_length
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        end_to_end_length_per_link = model.end_to_end_length_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            end_to_end_length_per_link / link_length -
            nondimensional_end_to_end_length_per_link
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
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
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        force = model.force(potential_distance, potential_stiffness)
        residual_abs =
            force / BOLTZMANN_CONSTANT / temperature * link_length - nondimensional_force
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        gibbs_free_energy =
            model.gibbs_free_energy(potential_distance, potential_stiffness, temperature)
        residual_abs =
            gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            relative_gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::nondimensional::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            relative_gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        end_to_end_length =
            model.end_to_end_length(potential_distance, potential_stiffness, temperature)
        end_to_end_length_per_link = model.end_to_end_length_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs = end_to_end_length / number_of_links - end_to_end_length_per_link
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
        )
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_end_to_end_length / number_of_links -
            nondimensional_end_to_end_length_per_link
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        gibbs_free_energy =
            model.gibbs_free_energy(potential_distance, potential_stiffness, temperature)
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs = gibbs_free_energy / number_of_links - gibbs_free_energy_per_link
        residual_rel = residual_abs / gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            relative_gibbs_free_energy / number_of_links -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        residual_abs =
            nondimensional_gibbs_free_energy / number_of_links -
            nondimensional_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::per_link::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_relative_gibbs_free_energy / number_of_links -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::relative::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        gibbs_free_energy =
            model.gibbs_free_energy(potential_distance, potential_stiffness, temperature)
        gibbs_free_energy_0 = model.gibbs_free_energy(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs = gibbs_free_energy - gibbs_free_energy_0 - relative_gibbs_free_energy
        residual_rel = residual_abs / gibbs_free_energy_0
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::relative::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        residual_abs =
            gibbs_free_energy_per_link - gibbs_free_energy_per_link_0 -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / gibbs_free_energy_per_link_0
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::relative::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(
            ZERO,
            nondimensional_potential_stiffness,
            temperature,
        )
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_gibbs_free_energy - nondimensional_gibbs_free_energy_0 -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_0
        @test abs(residual_abs) <= number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::relative::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_0 =
            model.nondimensional_gibbs_free_energy_per_link(
                ZERO,
                nondimensional_potential_stiffness,
                temperature,
            )
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        residual_abs =
            nondimensional_gibbs_free_energy_per_link -
            nondimensional_gibbs_free_energy_per_link_0 -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link_0
        @test abs(residual_abs) <= number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::force" begin
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
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        force_0 = model.force(ZERO, potential_stiffness)
        @test abs(force_0) <= potential_stiffness * ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::nondimensional_force" begin
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
        nondimensional_force_0 =
            model.nondimensional_force(ZERO, nondimensional_potential_stiffness)
        @test abs(nondimensional_force_0) <=
              number_of_links * nondimensional_potential_stiffness * ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::relative_gibbs_free_energy" begin
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
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::relative_gibbs_free_energy_per_link" begin
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
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(
            ZERO * number_of_links * link_length,
            potential_stiffness,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::nondimensional_relative_gibbs_free_energy" begin
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
        nondimensional_relative_gibbs_free_energy_0 =
            model.nondimensional_relative_gibbs_free_energy(
                ZERO,
                nondimensional_potential_stiffness,
            )
        @test abs(nondimensional_relative_gibbs_free_energy_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::zero::nondimensional_relative_gibbs_free_energy_per_link" begin
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
        nondimensional_relative_gibbs_free_energy_per_link_0 =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                ZERO,
                nondimensional_potential_stiffness,
            )
        @test abs(nondimensional_relative_gibbs_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::connection::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        end_to_end_length =
            model.end_to_end_length(potential_distance, potential_stiffness, temperature)
        h = parameters.rel_tol * number_of_links * link_length
        end_to_end_length_from_derivative =
            -1.0 / potential_stiffness * (
                model.relative_gibbs_free_energy(
                    potential_distance + 0.5 * h,
                    potential_stiffness,
                    temperature,
                ) - model.relative_gibbs_free_energy(
                    potential_distance - 0.5 * h,
                    potential_stiffness,
                    temperature,
                )
            ) / h
        residual_abs = end_to_end_length - end_to_end_length_from_derivative
        residual_rel = residual_abs / end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::connection::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance =
            nondimensional_potential_distance * number_of_links * link_length
        potential_stiffness =
            nondimensional_potential_stiffness / link_length^2 *
            BOLTZMANN_CONSTANT *
            temperature
        end_to_end_length_per_link = model.end_to_end_length_per_link(
            potential_distance,
            potential_stiffness,
            temperature,
        )
        h = parameters.rel_tol * number_of_links * link_length
        end_to_end_length_per_link_from_derivative =
            -1.0 / potential_stiffness * (
                model.relative_gibbs_free_energy_per_link(
                    potential_distance + 0.5 * h,
                    potential_stiffness,
                    temperature,
                ) - model.relative_gibbs_free_energy_per_link(
                    potential_distance - 0.5 * h,
                    potential_stiffness,
                    temperature,
                )
            ) / h
        residual_abs =
            end_to_end_length_per_link - end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::connection::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
        )
        h = parameters.rel_tol
        nondimensional_end_to_end_length_from_derivative =
            -1.0 / nondimensional_potential_stiffness / number_of_links * (
                model.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance + 0.5 * h,
                    nondimensional_potential_stiffness,
                ) - model.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance - 0.5 * h,
                    nondimensional_potential_stiffness,
                )
            ) / h
        residual_abs =
            nondimensional_end_to_end_length -
            nondimensional_end_to_end_length_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::test::connection::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_potential_distance =
            parameters.nondimensional_potential_distance_reference +
            parameters.nondimensional_potential_distance_scale * (0.5 - rand())
        nondimensional_potential_stiffness =
            parameters.nondimensional_potential_stiffness_reference +
            parameters.nondimensional_potential_stiffness_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            )
        h = parameters.rel_tol
        nondimensional_end_to_end_length_per_link_from_derivative =
            -1.0 / nondimensional_potential_stiffness / number_of_links * (
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_potential_distance + 0.5 * h,
                    nondimensional_potential_stiffness,
                ) - model.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_potential_distance - 0.5 * h,
                    nondimensional_potential_stiffness,
                )
            ) / h
        residual_abs =
            nondimensional_end_to_end_length_per_link -
            nondimensional_end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

end
