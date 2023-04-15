module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT, PLANCK_CONSTANT
using Polymers.Physics.SingleChain: ZERO, POINTS, integrate, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics: FJC

@testset "physics::single_chain::fjc::thermodynamics::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::test::base::number_of_links" begin
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

@testset "physics::single_chain::fjc::thermodynamics::test::base::link_length" begin
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

@testset "physics::single_chain::fjc::thermodynamics::test::base::hinge_mass" begin
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

@testset "physics::single_chain::fjc::thermodynamics::test::base::all_parameters" begin
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

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::force" begin
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        force_out = model.isometric.legendre.force(end_to_end_length, temperature)
        residual_abs = force - force_out
        residual_rel = residual_abs / force
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::nondimensional_force" begin
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
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_force_out = model.isometric.legendre.nondimensional_force(
            nondimensional_end_to_end_length_per_link,
        )
        residual_abs = nondimensional_force - nondimensional_force_out
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::helmholtz_free_energy" begin
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        helmholtz_free_energy_legendre =
            model.isotensional.gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        helmholtz_free_energy_legendre_out =
            model.isometric.legendre.helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy_legendre - helmholtz_free_energy_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            log(
                8.0 * pi^2 * hinge_mass * link_length^2 * BOLTZMANN_CONSTANT * temperature /
                PLANCK_CONSTANT^2,
            )
        residual_rel = residual_abs / helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::helmholtz_free_energy_per_link" begin
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.end_to_end_length_per_link(force, temperature)
        helmholtz_free_energy_per_link_legendre =
            model.isotensional.gibbs_free_energy_per_link(force, temperature) +
            force * end_to_end_length_per_link
        helmholtz_free_energy_per_link_legendre_out =
            model.isometric.legendre.helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            helmholtz_free_energy_per_link_legendre -
            helmholtz_free_energy_per_link_legendre_out +
            BOLTZMANN_CONSTANT *
            temperature *
            log(
                8.0 * pi^2 * hinge_mass * link_length^2 * BOLTZMANN_CONSTANT * temperature /
                PLANCK_CONSTANT^2,
            ) / number_of_links
        residual_rel = residual_abs / helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::relative_helmholtz_free_energy" begin
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        relative_helmholtz_free_energy_legendre =
            model.isotensional.relative_gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        relative_helmholtz_free_energy_legendre_out =
            model.isometric.legendre.relative_helmholtz_free_energy(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_legendre -
            relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::relative_helmholtz_free_energy_per_link" begin
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_per_link =
            model.isotensional.end_to_end_length_per_link(force, temperature)
        relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.relative_gibbs_free_energy_per_link(force, temperature) +
            force * end_to_end_length_per_link
        relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.legendre.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        residual_abs =
            relative_helmholtz_free_energy_per_link_legendre -
            relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::nondimensional_helmholtz_free_energy" begin
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
            model.isotensional.nondimensional_end_to_end_length(nondimensional_force)
        nondimensional_end_to_end_length_per_link =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_helmholtz_free_energy_legendre =
            model.isotensional.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_helmholtz_free_energy_legendre_out =
            model.isometric.legendre.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_legendre -
            nondimensional_helmholtz_free_energy_legendre_out + log(
                8.0 * pi^2 * hinge_mass * link_length^2 * BOLTZMANN_CONSTANT * temperature /
                PLANCK_CONSTANT^2,
            )
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::nondimensional_helmholtz_free_energy_per_link" begin
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
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_helmholtz_free_energy_per_link_legendre =
            model.isotensional.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.legendre.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link_legendre -
            nondimensional_helmholtz_free_energy_per_link_legendre_out + log(
                8.0 * pi^2 * hinge_mass * link_length^2 * BOLTZMANN_CONSTANT * temperature /
                PLANCK_CONSTANT^2,
            ) / number_of_links
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::nondimensional_relative_helmholtz_free_energy" begin
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
            model.isotensional.nondimensional_end_to_end_length(nondimensional_force)
        nondimensional_end_to_end_length_per_link =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_relative_helmholtz_free_energy_legendre =
            model.isotensional.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
            ) + nondimensional_force * nondimensional_end_to_end_length
        nondimensional_relative_helmholtz_free_energy_legendre_out =
            model.isometric.legendre.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_legendre -
            nondimensional_relative_helmholtz_free_energy_legendre_out
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::legendre::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_legendre =
            model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        nondimensional_relative_helmholtz_free_energy_per_link_legendre_out =
            model.isometric.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link_legendre -
            nondimensional_relative_helmholtz_free_energy_per_link_legendre_out
        residual_rel =
            residual_abs / nondimensional_relative_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        end_to_end_length_out = model.isotensional.end_to_end_length(force, temperature)
        residual_abs = end_to_end_length - end_to_end_length_out
        residual_rel = residual_abs / end_to_end_length
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        end_to_end_length_per_link = nondimensional_end_to_end_length_per_link * link_length
        end_to_end_length_per_link_out =
            model.isotensional.end_to_end_length_per_link(force, temperature)
        residual_abs = end_to_end_length_per_link - end_to_end_length_per_link_out
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links
        nondimensional_end_to_end_length_out =
            model.isotensional.nondimensional_end_to_end_length(nondimensional_force)
        residual_abs =
            nondimensional_end_to_end_length - nondimensional_end_to_end_length_out
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_end_to_end_length_per_link_out =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link -
            nondimensional_end_to_end_length_per_link_out
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        force_out = model.isometric.force(end_to_end_length, temperature)
        residual_abs = force - force_out
        residual_rel = residual_abs / force
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_force_out =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        residual_abs = nondimensional_force - nondimensional_force_out
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        helmholtz_free_energy =
            model.isometric.helmholtz_free_energy(end_to_end_length, temperature)
        helmholtz_free_energy_out =
            model.isotensional.gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        residual_abs = helmholtz_free_energy - helmholtz_free_energy_out
        residual_rel = residual_abs / helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        end_to_end_length_per_link = nondimensional_end_to_end_length_per_link * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        helmholtz_free_energy_per_link =
            model.isometric.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        helmholtz_free_energy_per_link_out =
            model.isotensional.gibbs_free_energy_per_link(force, temperature) +
            force * end_to_end_length_per_link
        residual_abs = helmholtz_free_energy_per_link - helmholtz_free_energy_per_link_out
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        relative_helmholtz_free_energy =
            model.isometric.relative_helmholtz_free_energy(end_to_end_length, temperature)
        relative_helmholtz_free_energy_out =
            model.isotensional.relative_gibbs_free_energy(force, temperature) +
            force * end_to_end_length
        residual_abs = relative_helmholtz_free_energy - relative_helmholtz_free_energy_out
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        end_to_end_length_per_link = nondimensional_end_to_end_length_per_link * link_length
        force = model.isometric.force(end_to_end_length, temperature)
        relative_helmholtz_free_energy_per_link =
            model.isometric.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            )
        relative_helmholtz_free_energy_per_link_out =
            model.isotensional.relative_gibbs_free_energy_per_link(force, temperature) +
            force * end_to_end_length_per_link
        residual_abs =
            relative_helmholtz_free_energy_per_link -
            relative_helmholtz_free_energy_per_link_out
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_helmholtz_free_energy =
            model.isometric.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_out =
            model.isotensional.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length
        residual_abs =
            nondimensional_helmholtz_free_energy - nondimensional_helmholtz_free_energy_out
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_helmholtz_free_energy_per_link =
            model.isometric.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_out =
            model.isotensional.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link -
            nondimensional_helmholtz_free_energy_per_link_out
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_relative_helmholtz_free_energy =
            model.isometric.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_out =
            model.isotensional.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
            ) + nondimensional_force * nondimensional_end_to_end_length
        residual_abs =
            nondimensional_relative_helmholtz_free_energy -
            nondimensional_relative_helmholtz_free_energy_out
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_force =
            model.isometric.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_out =
            model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
            ) + nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link -
            nondimensional_relative_helmholtz_free_energy_per_link_out
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        gibbs_free_energy = model.isotensional.gibbs_free_energy(force, temperature)
        gibbs_free_energy_out =
            model.isometric.helmholtz_free_energy(end_to_end_length, temperature) -
            force * end_to_end_length
        residual_abs = gibbs_free_energy - gibbs_free_energy_out
        residual_rel = residual_abs / gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_per_link = end_to_end_length / number_of_links
        gibbs_free_energy_per_link =
            model.isotensional.gibbs_free_energy_per_link(force, temperature)
        gibbs_free_energy_per_link_out =
            model.isometric.helmholtz_free_energy_per_link(end_to_end_length, temperature) -
            force * end_to_end_length_per_link
        residual_abs = gibbs_free_energy_per_link - gibbs_free_energy_per_link_out
        residual_rel = residual_abs / gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        relative_gibbs_free_energy =
            model.isotensional.relative_gibbs_free_energy(force, temperature)
        relative_gibbs_free_energy_out =
            model.isometric.relative_helmholtz_free_energy(end_to_end_length, temperature) -
            force * end_to_end_length
        residual_abs = relative_gibbs_free_energy - relative_gibbs_free_energy_out
        residual_rel = residual_abs / relative_gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
        end_to_end_length = model.isotensional.end_to_end_length(force, temperature)
        end_to_end_length_per_link = end_to_end_length / number_of_links
        relative_gibbs_free_energy_per_link =
            model.isotensional.relative_gibbs_free_energy_per_link(force, temperature)
        relative_gibbs_free_energy_per_link_out =
            model.isometric.relative_helmholtz_free_energy_per_link(
                end_to_end_length,
                temperature,
            ) - force * end_to_end_length_per_link
        residual_abs =
            relative_gibbs_free_energy_per_link - relative_gibbs_free_energy_per_link_out
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links
        nondimensional_gibbs_free_energy =
            model.isotensional.nondimensional_gibbs_free_energy(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_out =
            model.isometric.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            ) - nondimensional_force * nondimensional_end_to_end_length
        residual_abs =
            nondimensional_gibbs_free_energy - nondimensional_gibbs_free_energy_out
        residual_rel = residual_abs / nondimensional_gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
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
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_gibbs_free_energy_per_link =
            model.isotensional.nondimensional_gibbs_free_energy_per_link(
                nondimensional_force,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_out =
            model.isometric.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            ) - nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs =
            nondimensional_gibbs_free_energy_per_link -
            nondimensional_gibbs_free_energy_per_link_out
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links
        nondimensional_relative_gibbs_free_energy =
            model.isotensional.nondimensional_relative_gibbs_free_energy(
                nondimensional_force,
            )
        nondimensional_relative_gibbs_free_energy_out =
            model.isometric.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            ) - nondimensional_force * nondimensional_end_to_end_length
        residual_abs =
            nondimensional_relative_gibbs_free_energy -
            nondimensional_relative_gibbs_free_energy_out
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_force =
            parameters.nondimensional_force_reference +
            parameters.nondimensional_force_scale * (0.5 - rand())
        nondimensional_end_to_end_length_per_link =
            model.isotensional.nondimensional_end_to_end_length_per_link(
                nondimensional_force,
            )
        nondimensional_relative_gibbs_free_energy_per_link =
            model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_force,
            )
        nondimensional_relative_gibbs_free_energy_per_link_out =
            model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            ) - nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs =
            nondimensional_relative_gibbs_free_energy_per_link -
            nondimensional_relative_gibbs_free_energy_per_link_out
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol_thermodynamic_limit
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.modified_canonical.force(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    ) - model.isometric.force(end_to_end_length, temperature)
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.modified_canonical.force(
                    end_to_end_length,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.modified_canonical.nondimensional_force(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    ) - model.isometric.nondimensional_force(
                        nondimensional_end_to_end_length_per_link,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.modified_canonical.nondimensional_force(
                    nondimensional_end_to_end_length_per_link,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.modified_canonical.relative_helmholtz_free_energy(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    ) - model.isometric.relative_helmholtz_free_energy(
                        end_to_end_length,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.modified_canonical.relative_helmholtz_free_energy(
                    end_to_end_length,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.modified_canonical.relative_helmholtz_free_energy_per_link(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    ) - model.isometric.relative_helmholtz_free_energy_per_link(
                        end_to_end_length,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.modified_canonical.relative_helmholtz_free_energy_per_link(
                    end_to_end_length,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.modified_canonical.nondimensional_relative_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    ) - model.isometric.nondimensional_relative_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.modified_canonical.nondimensional_relative_helmholtz_free_energy(
                    nondimensional_end_to_end_length_per_link,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_strong_potential_isometric::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    ) -
                    model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.end_to_end_length(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.end_to_end_length(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.end_to_end_length(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.end_to_end_length_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.end_to_end_length_per_link(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.end_to_end_length_per_link(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_end_to_end_length(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_end_to_end_length(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.nondimensional_end_to_end_length(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_end_to_end_length_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_end_to_end_length_per_link(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.nondimensional_end_to_end_length_per_link(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1 *
            number_of_links *
            link_length
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            force_ref = model.modified_canonical.force(
                potential_distance_ref,
                potential_stiffness,
                temperature,
            )
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.gibbs_free_energy(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.gibbs_free_energy(force, temperature) +
                    model.isotensional.gibbs_free_energy(force_ref, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return (
                    model.modified_canonical.gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.gibbs_free_energy(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1 *
            number_of_links *
            link_length
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            force_ref = model.modified_canonical.force(
                potential_distance_ref,
                potential_stiffness,
                temperature,
            )
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.gibbs_free_energy_per_link(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.gibbs_free_energy_per_link(force, temperature) +
                    model.isotensional.gibbs_free_energy_per_link(force_ref, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return (
                    model.modified_canonical.gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.gibbs_free_energy_per_link(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.relative_gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.relative_gibbs_free_energy(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.relative_gibbs_free_energy(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.force(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )
                return (
                    model.modified_canonical.relative_gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.relative_gibbs_free_energy_per_link(
                        force,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.relative_gibbs_free_energy_per_link(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            nondimensional_force_ref = model.modified_canonical.nondimensional_force(
                nondimensional_potential_distance_ref,
                nondimensional_potential_stiffness,
            )
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.isotensional.nondimensional_gibbs_free_energy(
                        nondimensional_force,
                        temperature,
                    ) + model.isotensional.nondimensional_gibbs_free_energy(
                        nondimensional_force_ref,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return (
                    model.modified_canonical.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            nondimensional_force_ref = model.modified_canonical.nondimensional_force(
                nondimensional_potential_distance_ref,
                nondimensional_potential_stiffness,
            )
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.isotensional.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_force,
                        temperature,
                    ) + model.isotensional.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_force_ref,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return (
                    model.modified_canonical.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.modified_canonical.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_relative_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_relative_gibbs_free_energy(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_weak_potential_isotensional::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force = model.modified_canonical.nondimensional_force(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )
                return (
                    model.modified_canonical.nondimensional_relative_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_relative_gibbs_free_energy(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.isometric.force(end_to_end_length, temperature) -
                    model.modified_canonical.asymptotic.strong_potential.force(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.isometric.force(end_to_end_length, temperature)^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_force(
                        nondimensional_end_to_end_length_per_link,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_force(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.isometric.nondimensional_force(
                    nondimensional_end_to_end_length_per_link,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.isometric.helmholtz_free_energy(end_to_end_length, temperature) -
                    model.isometric.helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    ) +
                    model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return (
                    model.isometric.helmholtz_free_energy(end_to_end_length, temperature) -
                    model.isometric.helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.isometric.helmholtz_free_energy_per_link(
                        end_to_end_length,
                        temperature,
                    ) - model.isometric.helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy_per_link(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    ) +
                    model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return (
                    model.isometric.helmholtz_free_energy_per_link(
                        end_to_end_length,
                        temperature,
                    ) - model.isometric.helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small *
                        number_of_links *
                        link_length,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.isometric.relative_helmholtz_free_energy(
                        end_to_end_length,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.isometric.relative_helmholtz_free_energy(
                    end_to_end_length,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(end_to_end_length)
                return (
                    model.isometric.relative_helmholtz_free_energy_per_link(
                        end_to_end_length,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy_per_link(
                        end_to_end_length,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(end_to_end_length)
                return model.isometric.relative_helmholtz_free_energy_per_link(
                    end_to_end_length,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO * number_of_links * link_length,
                parameters.nondimensional_potential_distance_small *
                number_of_links *
                link_length,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                        temperature,
                    ) - model.isometric.nondimensional_helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) +
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                        temperature,
                    ) - model.isometric.nondimensional_helmholtz_free_energy(
                        parameters.nondimensional_potential_distance_small,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                        temperature,
                    ) - model.isometric.nondimensional_helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) +
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                        temperature,
                    ) - model.isometric.nondimensional_helmholtz_free_energy_per_link(
                        parameters.nondimensional_potential_distance_small,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_relative_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.isometric.nondimensional_relative_helmholtz_free_energy(
                    nondimensional_end_to_end_length_per_link,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_strong_potential_isometric::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_end_to_end_length_per_link)
                return (
                    model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                    ) -
                    model.modified_canonical.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy_per_link(
                        nondimensional_end_to_end_length_per_link,
                        nondimensional_potential_stiffness,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_end_to_end_length_per_link)
                return model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                ZERO,
                parameters.nondimensional_potential_distance_small,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_large * parameters.log_log_scale,
        )
        log_log_slope = log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.end_to_end_length(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.end_to_end_length(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.asymptotic.weak_potential.end_to_end_length(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.end_to_end_length_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.end_to_end_length_per_link(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.asymptotic.weak_potential.end_to_end_length_per_link(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_end_to_end_length(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_end_to_end_length_per_link(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1 *
            number_of_links *
            link_length
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            force_ref = model.modified_canonical.asymptotic.weak_potential.force(
                potential_distance_ref,
                potential_stiffness,
            )
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.gibbs_free_energy(force, temperature) +
                    model.isotensional.gibbs_free_energy(force_ref, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return (
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1 *
            number_of_links *
            link_length
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            force_ref = model.modified_canonical.asymptotic.weak_potential.force(
                potential_distance_ref,
                potential_stiffness,
            )
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.gibbs_free_energy_per_link(force, temperature) +
                    model.isotensional.gibbs_free_energy_per_link(force_ref, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return (
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(
                        potential_distance_ref,
                        potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.relative_gibbs_free_energy(force, temperature)
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            potential_stiffness =
                BOLTZMANN_CONSTANT * temperature / link_length^2 *
                nondimensional_potential_stiffness
            function integrand_numerator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                force = model.modified_canonical.asymptotic.weak_potential.force(
                    potential_distance,
                    potential_stiffness,
                )
                return (
                    model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(
                        potential_distance,
                        potential_stiffness,
                        temperature,
                    ) - model.isotensional.relative_gibbs_free_energy_per_link(
                        force,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                potential_distance =
                    number_of_links * link_length * nondimensional_potential_distance
                return model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(
                    potential_distance,
                    potential_stiffness,
                    temperature,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            nondimensional_force_ref =
                model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                    nondimensional_potential_distance_ref,
                    nondimensional_potential_stiffness,
                )
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.isotensional.nondimensional_gibbs_free_energy(
                        nondimensional_force,
                        temperature,
                    ) + model.isotensional.nondimensional_gibbs_free_energy(
                        nondimensional_force_ref,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_potential_distance_ref =
            parameters.nondimensional_potential_distance_large_1
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            nondimensional_force_ref =
                model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                    nondimensional_potential_distance_ref,
                    nondimensional_potential_stiffness,
                )
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) - model.isotensional.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_force,
                        temperature,
                    ) + model.isotensional.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_force_ref,
                        temperature,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                        temperature,
                    ) -
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(
                        nondimensional_potential_distance_ref,
                        nondimensional_potential_stiffness,
                        temperature,
                    )
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_relative_gibbs_free_energy(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::test::modified_canonical_asymptotic_weak_potential_isotensional::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            parameters.number_of_links_maximum - parameters.number_of_links_minimum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        function residual_rel(nondimensional_potential_stiffness)
            function integrand_numerator(nondimensional_potential_distance)
                nondimensional_force =
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_force(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    )
                return (
                    model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(
                        nondimensional_potential_distance,
                        nondimensional_potential_stiffness,
                    ) - model.isotensional.nondimensional_relative_gibbs_free_energy(
                        nondimensional_force,
                    )
                )^2
            end
            function integrand_denominator(nondimensional_potential_distance)
                return model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(
                    nondimensional_potential_distance,
                    nondimensional_potential_stiffness,
                )^2
            end
            numerator = integrate(
                integrand_numerator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            denominator = integrate(
                integrand_denominator,
                parameters.nondimensional_potential_distance_large_1,
                parameters.nondimensional_potential_distance_large_2,
                POINTS,
            )
            return sqrt(numerator / denominator)
        end
        residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small)
        residual_rel_2 = residual_rel(
            parameters.nondimensional_potential_stiffness_small * parameters.log_log_scale,
        )
        log_log_slope =
            -log(residual_rel_2 / residual_rel_1) / log(parameters.log_log_scale)
        @test abs(log_log_slope + 1.0) <= parameters.log_log_tol
    end
end

end
