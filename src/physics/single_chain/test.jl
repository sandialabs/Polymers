module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Ideal: IDEAL
using Polymers.Physics.SingleChain.Fjc: FJC
using Polymers.Physics.SingleChain.Efjc: EFJC
using Polymers.Physics.SingleChain.Swfjc: SWFJC

@testset "physics::single_chain::test::fjc_ideal::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_fjc =
            fjc.thermodynamics.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length_fjc - end_to_end_length_ideal
        residual_rel = residual_abs / end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link_fjc =
            fjc.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length_per_link_fjc - end_to_end_length_per_link_ideal
        residual_rel = residual_abs / end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_fjc =
            fjc.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_fjc - nondimensional_end_to_end_length_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_per_link_fjc =
            fjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link_fjc -
            nondimensional_end_to_end_length_per_link_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force_fjc = fjc.thermodynamics.isometric.force(end_to_end_length, temperature)
        force_ideal = ideal.thermodynamics.isometric.force(end_to_end_length, temperature)
        residual_abs = force_fjc - force_ideal
        residual_rel = residual_abs / force_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nodimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        nondimensional_force_fjc = fjc.thermodynamics.isometric.nondimensional_force(
            nondimensional_end_to_end_length_per_link,
        )
        nondimensional_force_ideal = ideal.thermodynamics.isometric.nondimensional_force(
            nondimensional_end_to_end_length_per_link,
        )
        residual_abs = nondimensional_force_fjc - nondimensional_force_ideal
        residual_rel = residual_abs / nondimensional_force_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy_fjc = fjc.thermodynamics.isometric.helmholtz_free_energy(
            end_to_end_length,
            temperature,
        )
        helmholtz_free_energy_ideal = ideal.thermodynamics.isometric.helmholtz_free_energy(
            end_to_end_length,
            temperature,
        )
        residual_abs = helmholtz_free_energy_fjc - helmholtz_free_energy_ideal
        residual_rel = residual_abs / helmholtz_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy_per_link_fjc =
            fjc.thermodynamics.isometric.helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        helmholtz_free_energy_per_link_ideal =
            ideal.thermodynamics.isometric.helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_per_link_fjc - helmholtz_free_energy_per_link_ideal
        residual_rel = residual_abs / helmholtz_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy_fjc =
            fjc.thermodynamics.isometric.relative_helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        relative_helmholtz_free_energy_ideal =
            ideal.thermodynamics.isometric.relative_helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_fjc - relative_helmholtz_free_energy_ideal
        residual_rel = residual_abs / relative_helmholtz_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy_per_link_fjc =
            fjc.thermodynamics.isometric.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        relative_helmholtz_free_energy_per_link_ideal =
            ideal.thermodynamics.isometric.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_per_link_fjc -
            relative_helmholtz_free_energy_per_link_ideal
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_fjc =
            fjc.thermodynamics.isometric.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_ideal =
            ideal.thermodynamics.isometric.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_fjc -
            nondimensional_helmholtz_free_energy_ideal
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_per_link_fjc =
            fjc.thermodynamics.isometric.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_ideal =
            ideal.thermodynamics.isometric.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link_fjc -
            nondimensional_helmholtz_free_energy_per_link_ideal
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        nondimensional_relative_helmholtz_free_energy_fjc =
            fjc.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_ideal =
            ideal.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_fjc -
            nondimensional_relative_helmholtz_free_energy_ideal
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        nondimensional_relative_helmholtz_free_energy_per_link_fjc =
            fjc.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_ideal =
            ideal.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link_fjc -
            nondimensional_relative_helmholtz_free_energy_per_link_ideal
        residual_rel =
            residual_abs / nondimensional_relative_helmholtz_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * number_of_links * link_length
        gibbs_free_energy_fjc =
            fjc.thermodynamics.isotensional.gibbs_free_energy(force, temperature)
        gibbs_free_energy_ideal =
            ideal.thermodynamics.isotensional.gibbs_free_energy(force, temperature)
        residual_abs = gibbs_free_energy_fjc - gibbs_free_energy_ideal
        residual_rel = residual_abs / gibbs_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * number_of_links * link_length
        gibbs_free_energy_per_link_fjc =
            fjc.thermodynamics.isotensional.gibbs_free_energy_per_link(force, temperature)
        gibbs_free_energy_per_link_ideal =
            ideal.thermodynamics.isotensional.gibbs_free_energy_per_link(force, temperature)
        residual_abs = gibbs_free_energy_per_link_fjc - gibbs_free_energy_per_link_ideal
        residual_rel = residual_abs / gibbs_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * number_of_links * link_length
        relative_gibbs_free_energy_fjc =
            fjc.thermodynamics.isotensional.relative_gibbs_free_energy(force, temperature)
        relative_gibbs_free_energy_ideal =
            ideal.thermodynamics.isotensional.relative_gibbs_free_energy(force, temperature)
        residual_abs = relative_gibbs_free_energy_fjc - relative_gibbs_free_energy_ideal
        residual_rel = residual_abs / relative_gibbs_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * number_of_links * link_length
        relative_gibbs_free_energy_per_link_fjc =
            fjc.thermodynamics.isotensional.relative_gibbs_free_energy_per_link(
                force,
                temperature,
            )
        relative_gibbs_free_energy_per_link_ideal =
            ideal.thermodynamics.isotensional.relative_gibbs_free_energy_per_link(
                force,
                temperature,
            )
        residual_abs =
            relative_gibbs_free_energy_per_link_fjc -
            relative_gibbs_free_energy_per_link_ideal
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_fjc =
            fjc.thermodynamics.isotensional.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_ideal =
            ideal.thermodynamics.isotensional.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            )
        residual_abs =
            nondimensional_gibbs_free_energy_fjc - nondimensional_gibbs_free_energy_ideal
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link_fjc =
            fjc.thermodynamics.isotensional.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_ideal =
            ideal.thermodynamics.isotensional.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        residual_abs =
            nondimensional_gibbs_free_energy_per_link_fjc -
            nondimensional_gibbs_free_energy_per_link_ideal
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_relative_gibbs_free_energy_fjc =
            fjc.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
            )
        nondimensional_relative_gibbs_free_energy_ideal =
            ideal.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_relative_gibbs_free_energy_fjc -
            nondimensional_relative_gibbs_free_energy_ideal
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_relative_gibbs_free_energy_per_link_fjc =
            fjc.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
            )
        nondimensional_relative_gibbs_free_energy_per_link_ideal =
            ideal.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_relative_gibbs_free_energy_per_link_fjc -
            nondimensional_relative_gibbs_free_energy_per_link_ideal
        residual_rel =
            residual_abs / nondimensional_relative_gibbs_free_energy_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::fjc_ideal::equilibrium_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        equilibrium_distribution_fjc =
            fjc.thermodynamics.isometric.equilibrium_distribution(end_to_end_length)
        equilibrium_distribution_ideal =
            ideal.thermodynamics.isometric.equilibrium_distribution(end_to_end_length)
        residual_abs = equilibrium_distribution_fjc - equilibrium_distribution_ideal
        residual_rel = residual_abs / equilibrium_distribution_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_equilibrium_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        nondimensional_equilibrium_distribution_fjc =
            fjc.thermodynamics.isometric.nondimensional_equilibrium_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_equilibrium_distribution_ideal =
            ideal.thermodynamics.isometric.nondimensional_equilibrium_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_equilibrium_distribution_fjc -
            nondimensional_equilibrium_distribution_ideal
        residual_rel = residual_abs / nondimensional_equilibrium_distribution_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        equilibrium_radial_distribution_fjc =
            fjc.thermodynamics.isometric.equilibrium_radial_distribution(end_to_end_length)
        equilibrium_radial_distribution_ideal =
            ideal.thermodynamics.isometric.equilibrium_radial_distribution(
                end_to_end_length,
            )
        residual_abs =
            equilibrium_radial_distribution_fjc - equilibrium_radial_distribution_ideal
        residual_rel = residual_abs / equilibrium_radial_distribution_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::fjc_ideal::nondimensional_equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        fjc = FJC(number_of_links, link_length, hinge_mass)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_small *
            (1.0 - 0.5 * rand())
        nondimensional_equilibrium_radial_distribution_fjc =
            fjc.thermodynamics.isometric.nondimensional_equilibrium_radial_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_equilibrium_radial_distribution_ideal =
            ideal.thermodynamics.isometric.nondimensional_equilibrium_radial_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_equilibrium_radial_distribution_fjc -
            nondimensional_equilibrium_radial_distribution_ideal
        residual_rel = residual_abs / nondimensional_equilibrium_radial_distribution_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_end_to_end_length_per_link
    end
end

@testset "physics::single_chain::test::efjc_ideal::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        link_stiffness = parameters.link_stiffness_scale
        efjc = EFJC(number_of_links, link_length, hinge_mass, link_stiffness)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_efjc =
            efjc.thermodynamics.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length_efjc - end_to_end_length_ideal
        residual_rel = residual_abs / end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::efjc_ideal::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        link_stiffness = parameters.link_stiffness_scale
        efjc = EFJC(number_of_links, link_length, hinge_mass, link_stiffness)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link_efjc =
            efjc.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length_per_link_efjc - end_to_end_length_per_link_ideal
        residual_rel = residual_abs / end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::efjc_ideal::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        link_stiffness = parameters.link_stiffness_scale
        efjc = EFJC(number_of_links, link_length, hinge_mass, link_stiffness)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length_efjc =
            efjc.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_efjc - nondimensional_end_to_end_length_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::efjc_ideal::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        link_stiffness = parameters.link_stiffness_scale
        efjc = EFJC(number_of_links, link_length, hinge_mass, link_stiffness)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link_efjc =
            efjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link_efjc -
            nondimensional_end_to_end_length_per_link_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::swfjc_ideal::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_swfjc =
            swfjc.thermodynamics.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length_swfjc - end_to_end_length_ideal
        residual_rel = residual_abs / end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::swfjc_ideal::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link_swfjc =
            swfjc.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length_per_link_swfjc - end_to_end_length_per_link_ideal
        residual_rel = residual_abs / end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::swfjc_ideal::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_swfjc =
            swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_swfjc - nondimensional_end_to_end_length_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::swfjc_ideal::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        ideal = IDEAL(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_per_link_swfjc =
            swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_per_link_ideal =
            ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link_swfjc -
            nondimensional_end_to_end_length_per_link_ideal
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link_ideal
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit &&
              abs(residual_rel) <= nondimensional_force
    end
end

@testset "physics::single_chain::test::swfjc_fjc::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        fjc = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_swfjc =
            swfjc.thermodynamics.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_fjc =
            fjc.thermodynamics.isotensional.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length_swfjc - end_to_end_length_fjc
        @test abs(residual_abs) <= number_of_links * well_width
    end
end

@testset "physics::single_chain::test::swfjc_fjc::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        fjc = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force = nondimensional_force * BOLTZMANN_CONSTANT * temperature / link_length
        end_to_end_length_per_link_swfjc =
            swfjc.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        end_to_end_length_per_link_fjc =
            fjc.thermodynamics.isotensional.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length_per_link_swfjc - end_to_end_length_per_link_fjc
        @test abs(residual_abs) <= well_width
    end
end

@testset "physics::single_chain::test::swfjc_fjc::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        fjc = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_swfjc =
            swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_fjc =
            fjc.thermodynamics.isotensional.nondimensional_end_to_end_length(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_swfjc - nondimensional_end_to_end_length_fjc
        @test abs(residual_abs) <=
              number_of_links * parameters.nondimensional_well_width_small
    end
end

@testset "physics::single_chain::test::swfjc_fjc::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width = parameters.nondimensional_well_width_small * link_length
        swfjc = SWFJC(number_of_links, link_length, hinge_mass, well_width)
        fjc = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force = parameters.nondimensional_force_small * (1.0 - 0.5 * rand())
        nondimensional_end_to_end_length_per_link_swfjc =
            swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_end_to_end_length_per_link_fjc =
            fjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link_swfjc -
            nondimensional_end_to_end_length_per_link_fjc
        @test abs(residual_abs) <= parameters.nondimensional_well_width_small
    end
end

end
