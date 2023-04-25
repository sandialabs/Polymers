module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT, PLANCK_CONSTANT
using Polymers.Physics.SingleChain: parameters
using Polymers.Physics.SingleChain.Swfjc.Thermodynamics: SWFJC

@testset "physics::single_chain::swfjc::thermodynamics::test::base::init" begin
    @test isa(
        SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::swfjc::thermodynamics::test::base::number_of_links" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        @test SWFJC(
            number_of_links,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ).number_of_links == number_of_links
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::base::link_length" begin
    for _ = 1:parameters.number_of_loops
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
            parameters.well_width_reference,
        ).link_length == link_length
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::base::hinge_mass" begin
    for _ = 1:parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
            parameters.well_width_reference,
        ).hinge_mass == hinge_mass
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::base::well_width" begin
    for _ = 1:parameters.number_of_loops
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        @test SWFJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
            well_width,
        ).well_width == well_width
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::base::all_parameters" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        @test all(
            SWFJC(number_of_links, link_length, hinge_mass, well_width).number_of_links ==
            number_of_links &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).link_length ==
            link_length &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).hinge_mass ==
            hinge_mass &&
            SWFJC(number_of_links, link_length, hinge_mass, well_width).well_width ==
            well_width,
        )
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= 3e1 * parameters.abs_tol || abs(residual_rel) <= 3e1 * parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= 3e1 * parameters.abs_tol || abs(residual_rel) <= 3e1 * parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
            nondimensional_helmholtz_free_energy_per_link_legendre_out +
            log(
                8.0 * pi^2 * hinge_mass * link_length^2 * BOLTZMANN_CONSTANT * temperature /
                PLANCK_CONSTANT^2,
            ) / number_of_links
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link_legendre
        @test abs(residual_abs) <= parameters.abs_tol ||
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= 3e1 * parameters.abs_tol || abs(residual_rel) <= 3e1 * parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::test::legendre::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        well_width =
            parameters.well_width_reference + parameters.well_width_scale * (0.5 - rand())
        model = SWFJC(number_of_links, link_length, hinge_mass, well_width)
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
        @test abs(residual_abs) <= 3e1 * parameters.abs_tol || abs(residual_rel) <= 3e1 * parameters.rel_tol
    end
end

end
