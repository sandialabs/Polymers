module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT, PLANCK_CONSTANT
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Wlc.Thermodynamics: WLC

@testset "physics::single_chain::wlc::thermodynamics::test::base::init" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::base::number_of_links" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::base::link_length" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::base::hinge_mass" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::base::persistance_length" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::base::all_parameters" begin
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_gibbs_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

@testset "physics::single_chain::wlc::thermodynamics::test::thermodynamic_limit::nondimensional_relative_gibbs_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links = parameters.number_of_links_maximum
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        persistance_length =
            parameters.nondimensional_persistance_length_small *
            number_of_links *
            link_length
        model = WLC(number_of_links, link_length, hinge_mass, persistance_length)
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

end
