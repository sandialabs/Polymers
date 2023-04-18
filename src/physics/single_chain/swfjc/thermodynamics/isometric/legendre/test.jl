module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ONE, ZERO, POINTS, integrate, parameters
using Polymers.Physics.SingleChain.Swfjc.Thermodynamics.Isometric.Legendre: SWFJC

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::init" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::number_of_links" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::link_length" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::hinge_mass" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::well_width" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::base::all_parameters" begin
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::normalization::equilibrium_distribution" begin
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
        normalization = integrate(
            end_to_end_length ->
                4.0 *
                pi *
                end_to_end_length^2 *
                model.equilibrium_distribution(end_to_end_length),
            ZERO,
            ONE * number_of_links * (link_length + well_width),
            POINTS,
        )
        @test abs(normalization - 1.0) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::normalization::equilibrium_radial_distribution" begin
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
        normalization = integrate(
            end_to_end_length ->
                model.equilibrium_radial_distribution(end_to_end_length),
            ZERO,
            ONE * number_of_links * (link_length + well_width),
            POINTS,
        )
        @test abs(normalization - 1.0) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::normalization::nondimensional_equilibrium_distribution" begin
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
        normalization = integrate(
            nondimensional_end_to_end_length_per_link ->
                4.0 *
                pi *
                nondimensional_end_to_end_length_per_link^2 *
                model.nondimensional_equilibrium_distribution(
                    nondimensional_end_to_end_length_per_link,
                ),
            ZERO,
            ONE * (1.0 + well_width / link_length),
            POINTS,
        )
        @test abs(normalization - 1.0) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::normalization::nondimensional_equilibrium_radial_distribution" begin
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
        normalization = integrate(
            nondimensional_end_to_end_length_per_link ->
                model.nondimensional_equilibrium_radial_distribution(
                    nondimensional_end_to_end_length_per_link,
                ),
            ZERO,
            ONE * (1.0 + well_width / link_length),
            POINTS,
        )
        @test abs(normalization - 1.0) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::nondimensional::force" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.force(end_to_end_length, temperature)
        residual_abs =
            force / BOLTZMANN_CONSTANT / temperature * link_length - nondimensional_force
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::nondimensional::helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy = model.helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::nondimensional::helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::nondimensional::relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs =
            relative_helmholtz_free_energy / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::nondimensional::relative_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            relative_helmholtz_free_energy_per_link / BOLTZMANN_CONSTANT / temperature -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::per_link::helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy = model.helmholtz_free_energy(end_to_end_length, temperature)
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy / number_of_links - helmholtz_free_energy_per_link
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::per_link::relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            relative_helmholtz_free_energy / number_of_links -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::per_link::nondimensional_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
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

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::per_link::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy / number_of_links -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::relative::helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy = model.helmholtz_free_energy(end_to_end_length, temperature)
        helmholtz_free_energy_0 =
            model.helmholtz_free_energy(ZERO * number_of_links * link_length, temperature)
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy - helmholtz_free_energy_0 - relative_helmholtz_free_energy
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::relative::helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(
            ZERO * number_of_links * link_length,
            temperature,
        )
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs =
            helmholtz_free_energy_per_link - helmholtz_free_energy_per_link_0 -
            relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::relative::nondimensional_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        nondimensional_helmholtz_free_energy_0 =
            model.nondimensional_helmholtz_free_energy(ZERO, temperature)
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy - nondimensional_helmholtz_free_energy_0 -
            nondimensional_relative_helmholtz_free_energy
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::relative::nondimensional_helmholtz_free_energy_per_link" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_0 =
            model.nondimensional_helmholtz_free_energy_per_link(ZERO, temperature)
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link -
            nondimensional_helmholtz_free_energy_per_link_0 -
            nondimensional_relative_helmholtz_free_energy_per_link
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::force" begin
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
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        force_0 = model.force(ZERO * number_of_links * link_length, temperature)
        @test abs(force_0) <=
              3.1 * ZERO * number_of_links * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::nondimensional_force" begin
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
        nondimensional_force_0 = model.nondimensional_force(ZERO)
        @test abs(nondimensional_force_0) <= 3.1 * ZERO * number_of_links
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::relative_helmholtz_free_energy" begin
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
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(
            ZERO * number_of_links * link_length,
            temperature,
        )
        @test abs(relative_helmholtz_free_energy_0) <=
              ZERO * number_of_links * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::relative_helmholtz_free_energy_per_link" begin
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
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        relative_helmholtz_free_energy_per_link_0 =
            model.relative_helmholtz_free_energy_per_link(
                ZERO * number_of_links * link_length,
                temperature,
            )
        @test abs(relative_helmholtz_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_relative_helmholtz_free_energy_0 =
            model.nondimensional_relative_helmholtz_free_energy(ZERO)
        @test abs(nondimensional_relative_helmholtz_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
        nondimensional_relative_helmholtz_free_energy_per_link_0 =
            model.nondimensional_relative_helmholtz_free_energy_per_link(ZERO)
        @test abs(nondimensional_relative_helmholtz_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::equilibrium_radial_distribution" begin
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
        equilibrium_radial_distribution_0 =
            model.equilibrium_radial_distribution(ZERO * number_of_links * link_length)
        @test abs(equilibrium_radial_distribution_0) <= ZERO
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::zero::nondimensional_equilibrium_radial_distribution" begin
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
        nondimensional_equilibrium_radial_distribution_0 =
            model.equilibrium_radial_distribution(ZERO)
        @test abs(nondimensional_equilibrium_radial_distribution_0) <= ZERO
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::connection::force" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        force = model.force(end_to_end_length, temperature)
        h = parameters.rel_tol * number_of_links * link_length
        force_from_derivative =
            (
                model.relative_helmholtz_free_energy(
                    end_to_end_length + 0.5 * h,
                    temperature,
                ) - model.relative_helmholtz_free_energy(
                    end_to_end_length - 0.5 * h,
                    temperature,
                )
            ) / h
        residual_abs = force - force_from_derivative
        residual_rel = residual_abs / force
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::connection::nondimensional_force" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        nondimensional_force =
            model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        h = parameters.rel_tol
        nondimensional_force_from_derivative =
            (
                model.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link + 0.5 * h,
                ) - model.nondimensional_relative_helmholtz_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link - 0.5 * h,
                )
            ) / h
        residual_abs = nondimensional_force - nondimensional_force_from_derivative
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::connection::relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference + parameters.temperature_scale * (0.5 - rand())
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        relative_helmholtz_free_energy_from_connection =
            BOLTZMANN_CONSTANT *
            temperature *
            log((
                model.equilibrium_distribution(ZERO * number_of_links * link_length) /
                model.equilibrium_distribution(end_to_end_length)
            ))
        residual_abs =
            relative_helmholtz_free_energy - relative_helmholtz_free_energy_from_connection
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::swfjc::thermodynamics::isometric::legendre::test::connection::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_from_connection = log((
            model.nondimensional_equilibrium_distribution(ZERO) /
            model.nondimensional_equilibrium_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        ))
        residual_abs =
            nondimensional_relative_helmholtz_free_energy -
            nondimensional_relative_helmholtz_free_energy_from_connection
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= parameters.rel_tol
    end
end

end
