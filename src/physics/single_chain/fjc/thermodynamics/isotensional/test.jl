module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ZERO, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isotensional: FJC

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::base::number_of_links" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::base::link_length" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::base::hinge_mass" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::base::all_parameters" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            model.nondimensional_end_to_end_length(nondimensional_force)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length = model.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length / link_length - nondimensional_end_to_end_length
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(nondimensional_force)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link = model.end_to_end_length_per_link(force, temperature)
        residual_abs =
            end_to_end_length_per_link / link_length -
            nondimensional_end_to_end_length_per_link
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy =
            model.nondimensional_gibbs_free_energy(nondimensional_force, temperature)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        gibbs_free_energy = model.gibbs_free_energy(force, temperature)
        residual_abs =
            gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(force, temperature)
        residual_abs =
            gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(nondimensional_force)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(force, temperature)
        residual_abs =
            relative_gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::nondimensional::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(force, temperature)
        residual_abs =
            relative_gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length = model.end_to_end_length(force, temperature)
        end_to_end_length_per_link = model.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length / number_of_links - end_to_end_length_per_link
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            model.nondimensional_end_to_end_length(nondimensional_force)
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(nondimensional_force)
        residual_abs =
            nondimensional_end_to_end_length / number_of_links -
            nondimensional_end_to_end_length_per_link
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        gibbs_free_energy = model.gibbs_free_energy(force, temperature)
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(force, temperature)
        residual_abs = gibbs_free_energy / number_of_links - gibbs_free_energy_per_link
        residual_rel = residual_abs / gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(force, temperature)
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(force, temperature)
        residual_abs =
            relative_gibbs_free_energy / number_of_links -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy =
            model.nondimensional_gibbs_free_energy(nondimensional_force, temperature)
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
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

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::per_link::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(nondimensional_force)
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
        residual_abs =
            nondimensional_relative_gibbs_free_energy / number_of_links -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::relative::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        gibbs_free_energy = model.gibbs_free_energy(force, temperature)
        gibbs_free_energy_0 = model.gibbs_free_energy(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(force, temperature)
        residual_abs = gibbs_free_energy - gibbs_free_energy_0 - relative_gibbs_free_energy
        residual_rel = residual_abs / gibbs_free_energy_0
        @test abs(residual_abs) <=
              BOLTZMANN_CONSTANT * temperature * number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::relative::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(force, temperature)
        gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(force, temperature)
        residual_abs =
            gibbs_free_energy_per_link - gibbs_free_energy_per_link_0 -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / gibbs_free_energy_per_link_0
        @test abs(residual_abs) <= BOLTZMANN_CONSTANT * temperature * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::relative::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy =
            model.nondimensional_gibbs_free_energy(nondimensional_force, temperature)
        nondimensional_gibbs_free_energy_0 =
            model.nondimensional_gibbs_free_energy(ZERO, temperature)
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(nondimensional_force)
        residual_abs =
            nondimensional_gibbs_free_energy - nondimensional_gibbs_free_energy_0 -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_0
        @test abs(residual_abs) <= number_of_links * parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::relative::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_0 =
            model.nondimensional_gibbs_free_energy_per_link(ZERO, temperature)
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
        residual_abs =
            nondimensional_gibbs_free_energy_per_link -
            nondimensional_gibbs_free_energy_per_link_0 -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link_0
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length_0 = model.end_to_end_length(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        @test abs(end_to_end_length_0) <= ZERO * number_of_links * link_length
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length_per_link_0 = model.end_to_end_length_per_link(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        @test abs(end_to_end_length_per_link_0) <= ZERO * link_length
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_0 = model.nondimensional_end_to_end_length(ZERO)
        @test abs(nondimensional_end_to_end_length_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link_0 =
            model.nondimensional_end_to_end_length_per_link(ZERO)
        @test abs(nondimensional_end_to_end_length_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature * number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_relative_gibbs_free_energy_0 =
            model.nondimensional_relative_gibbs_free_energy(ZERO)
        @test abs(nondimensional_relative_gibbs_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::zero::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_relative_gibbs_free_energy_per_link_0 =
            model.nondimensional_relative_gibbs_free_energy_per_link(ZERO)
        @test abs(nondimensional_relative_gibbs_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::connection::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length = model.end_to_end_length(force, temperature)
        h = parameters.rel_tol * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_from_derivative =
            -(
                model.relative_gibbs_free_energy(force + 0.5 * h, temperature) -
                model.relative_gibbs_free_energy(force - 0.5 * h, temperature)
            ) / h
        residual_abs = end_to_end_length - end_to_end_length_from_derivative
        residual_rel = residual_abs / end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::connection::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link = model.end_to_end_length_per_link(force, temperature)
        h = parameters.rel_tol * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link_from_derivative =
            -(
                model.relative_gibbs_free_energy_per_link(force + 0.5 * h, temperature) -
                model.relative_gibbs_free_energy_per_link(force - 0.5 * h, temperature)
            ) / h
        residual_abs =
            end_to_end_length_per_link - end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::connection::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            model.nondimensional_end_to_end_length(nondimensional_force)
        h = parameters.rel_tol
        nondimensional_end_to_end_length_from_derivative =
            -(
                model.nondimensional_relative_gibbs_free_energy(
                    nondimensional_force + 0.5 * h,
                ) - model.nondimensional_relative_gibbs_free_energy(
                    nondimensional_force - 0.5 * h,
                )
            ) / h
        residual_abs =
            nondimensional_end_to_end_length -
            nondimensional_end_to_end_length_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::connection::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(nondimensional_force)
        h = parameters.rel_tol
        nondimensional_end_to_end_length_per_link_from_derivative =
            -(
                model.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_force + 0.5 * h,
                ) - model.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_force - 0.5 * h,
                )
            ) / h
        residual_abs =
            nondimensional_end_to_end_length_per_link -
            nondimensional_end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length = model.end_to_end_length(force, temperature)
        gibbs_free_energy = model.gibbs_free_energy(force, temperature)
        gibbs_free_energy_legendre =
            model.legendre.helmholtz_free_energy(force, temperature) -
            force * end_to_end_length
        residual_abs = gibbs_free_energy - gibbs_free_energy_legendre
        residual_rel = residual_abs / gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link = model.end_to_end_length_per_link(force, temperature)
        gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(force, temperature)
        gibbs_free_energy_per_link_legendre =
            model.legendre.helmholtz_free_energy_per_link(force, temperature) -
            force * end_to_end_length_per_link
        residual_abs = gibbs_free_energy_per_link - gibbs_free_energy_per_link_legendre
        residual_rel = residual_abs / gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length = model.end_to_end_length(force, temperature)
        end_to_end_length_0 = model.end_to_end_length(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_gibbs_free_energy = model.relative_gibbs_free_energy(force, temperature)
        relative_gibbs_free_energy_legendre =
            model.legendre.relative_helmholtz_free_energy(force, temperature) -
            force * end_to_end_length +
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length * end_to_end_length_0
        residual_abs = relative_gibbs_free_energy - relative_gibbs_free_energy_legendre
        residual_rel = residual_abs / relative_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link = model.end_to_end_length_per_link(force, temperature)
        end_to_end_length_per_link_0 = model.end_to_end_length_per_link(
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length,
            temperature,
        )
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(force, temperature)
        relative_gibbs_free_energy_per_link_legendre =
            model.legendre.relative_helmholtz_free_energy_per_link(force, temperature) -
            force * end_to_end_length_per_link +
            ZERO * BOLTZMANN_CONSTANT * temperature / link_length *
            end_to_end_length_per_link_0
        residual_abs =
            relative_gibbs_free_energy_per_link -
            relative_gibbs_free_energy_per_link_legendre
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            model.nondimensional_end_to_end_length(nondimensional_force)
        nondimensional_gibbs_free_energy =
            model.nondimensional_gibbs_free_energy(nondimensional_force, temperature)
        nondimensional_gibbs_free_energy_legendre =
            model.legendre.nondimensional_helmholtz_free_energy(
                nondimensional_force,
                temperature,
            ) - nondimensional_force * nondimensional_end_to_end_length
        residual_abs =
            nondimensional_gibbs_free_energy - nondimensional_gibbs_free_energy_legendre
        residual_rel = residual_abs / nondimensional_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(nondimensional_force)
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_legendre =
            model.legendre.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) - nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs =
            nondimensional_gibbs_free_energy_per_link -
            nondimensional_gibbs_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            model.nondimensional_end_to_end_length(nondimensional_force)
        nondimensional_end_to_end_length_0 = model.nondimensional_end_to_end_length(ZERO)
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(nondimensional_force)
        nondimensional_relative_gibbs_free_energy_legendre =
            model.legendre.nondimensional_relative_helmholtz_free_energy(
                nondimensional_force,
            ) - nondimensional_force * nondimensional_end_to_end_length +
            ZERO * nondimensional_end_to_end_length_0
        residual_abs =
            nondimensional_relative_gibbs_free_energy -
            nondimensional_relative_gibbs_free_energy_legendre
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.nondimensional_end_to_end_length_per_link(nondimensional_force)
        nondimensional_end_to_end_length_per_link_0 =
            model.nondimensional_end_to_end_length_per_link(ZERO)
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
        nondimensional_relative_gibbs_free_energy_per_link_legendre =
            model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_force,
            ) - nondimensional_force * nondimensional_end_to_end_length_per_link +
            ZERO * nondimensional_end_to_end_length_per_link_0
        residual_abs =
            nondimensional_relative_gibbs_free_energy_per_link -
            nondimensional_relative_gibbs_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre_connection::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        h = parameters.rel_tol * BOLTZMANN_CONSTANT * temperature / link_length
        force_from_derivative =
            (
                model.legendre.relative_helmholtz_free_energy(
                    force + 0.5 * h,
                    temperature,
                ) -
                model.legendre.relative_helmholtz_free_energy(force - 0.5 * h, temperature)
            ) / (
                model.end_to_end_length(force + 0.5 * h, temperature) -
                model.end_to_end_length(force - 0.5 * h, temperature)
            )
        residual_abs = force - force_from_derivative
        residual_rel = residual_abs / force
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isotensional::test::legendre_connection::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        h = parameters.rel_tol
        nondimensional_force_from_derivative =
            (
                model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_force + 0.5 * h,
                ) - model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_force - 0.5 * h,
                )
            ) / (
                model.nondimensional_end_to_end_length_per_link(
                    nondimensional_force + 0.5 * h,
                ) - model.nondimensional_end_to_end_length_per_link(
                    nondimensional_force - 0.5 * h,
                )
            )
        residual_abs = nondimensional_force - nondimensional_force_from_derivative
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= h
    end
end

end
