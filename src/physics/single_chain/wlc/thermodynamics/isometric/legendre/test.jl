module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ZERO, parameters
using Polymers.Physics.SingleChain.Wlc.Thermodynamics.Isometric.Legendre: WLC

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::init" begin
    @test isa(
        WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test WLC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.persistance_length_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.persistance_length_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::persistance_length" begin
    for _ = 1:parameters.number_of_loops
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        @test WLC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            persistance_length,
        ).persistance_length == persistance_length
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        @test all(
            WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length,
            ).number_of_links == number_of_links &&
            WLC(number_of_links, link_length, hinge_mass, persistance_length).link_length ==
            link_length &&
            WLC(number_of_links, link_length, hinge_mass, persistance_length).hinge_mass ==
            hinge_mass &&
            WLC(
                number_of_links,
                link_length,
                hinge_mass,
                persistance_length,
            ).persistance_length == persistance_length,
        )
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::nondimensional::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        gibbs_free_energy = model.gibbs_free_energy(end_to_end_length, temperature)
        residual_abs =
            gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::nondimensional::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        gibbs_free_energy_per_link =
            model.gibbs_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::nondimensional::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_gibbs_free_energy =
            model.relative_gibbs_free_energy(end_to_end_length, temperature)
        residual_abs =
            relative_gibbs_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::nondimensional::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            relative_gibbs_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::per_link::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        gibbs_free_energy = model.gibbs_free_energy(end_to_end_length, temperature)
        gibbs_free_energy_per_link =
            model.gibbs_free_energy_per_link(end_to_end_length, temperature)
        residual_abs = gibbs_free_energy / number_of_links - gibbs_free_energy_per_link
        residual_rel = residual_abs / gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::per_link::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_gibbs_free_energy =
            model.relative_gibbs_free_energy(end_to_end_length, temperature)
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            relative_gibbs_free_energy / number_of_links -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::per_link::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
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

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::per_link::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_gibbs_free_energy / number_of_links -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::relative::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        gibbs_free_energy = model.gibbs_free_energy(end_to_end_length, temperature)
        gibbs_free_energy_0 =
            model.gibbs_free_energy(ZERO * number_of_links * link_length, temperature)
        relative_gibbs_free_energy =
            model.relative_gibbs_free_energy(end_to_end_length, temperature)
        residual_abs = gibbs_free_energy - gibbs_free_energy_0 - relative_gibbs_free_energy
        residual_rel = residual_abs / relative_gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::relative::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        gibbs_free_energy_per_link =
            model.gibbs_free_energy_per_link(end_to_end_length, temperature)
        gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(
            ZERO * number_of_links * link_length,
            temperature,
        )
        relative_gibbs_free_energy_per_link =
            model.relative_gibbs_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            gibbs_free_energy_per_link - gibbs_free_energy_per_link_0 -
            relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / relative_gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::relative::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        nondimensional_gibbs_free_energy_0 =
            model.nondimensional_gibbs_free_energy(ZERO, temperature)
        nondimensional_relative_gibbs_free_energy =
            model.nondimensional_relative_gibbs_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_gibbs_free_energy - nondimensional_gibbs_free_energy_0 -
            nondimensional_relative_gibbs_free_energy
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::relative::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_gibbs_free_energy_per_link =
            model.nondimensional_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_gibbs_free_energy_per_link_0 =
            model.nondimensional_gibbs_free_energy_per_link(ZERO, temperature)
        nondimensional_relative_gibbs_free_energy_per_link =
            model.nondimensional_relative_gibbs_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_gibbs_free_energy_per_link -
            nondimensional_gibbs_free_energy_per_link_0 -
            nondimensional_relative_gibbs_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_gibbs_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::zero::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(
            ZERO * number_of_links * link_length,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_0) <=
              ZERO * number_of_links * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::zero::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(
            ZERO * number_of_links * link_length,
            temperature,
        )
        @test abs(relative_gibbs_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::zero::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_relative_gibbs_free_energy_0 =
            model.nondimensional_relative_gibbs_free_energy(ZERO)
        @test abs(nondimensional_relative_gibbs_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::wlc::thermodynamics::isometric::legendre::test::zero::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.persistance_length_reference +
            parameters.persistance_length_scale * (0.5 - rand())
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
        nondimensional_relative_gibbs_free_energy_per_link_0 =
            model.nondimensional_relative_gibbs_free_energy_per_link(ZERO)
        @test abs(nondimensional_relative_gibbs_free_energy_per_link_0) <= ZERO
    end
end

end
