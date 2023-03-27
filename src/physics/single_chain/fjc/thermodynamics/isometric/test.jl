module Test

using Test
using Polymers.Physics: BOLTZMANN_CONSTANT
using Polymers.Physics.SingleChain: ONE, ZERO, POINTS, integrate, parameters
using Polymers.Physics.SingleChain.Fjc.Thermodynamics.Isometric: FJC

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::init" begin
    @test isa(
        FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            parameters.hinge_mass_reference,
        ),
        Any,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::number_of_links" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::link_length" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::hinge_mass" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::all_parameters" begin
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::equilibrium_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        normalization = integrate(
            end_to_end_length ->
                4.0 *
                pi *
                end_to_end_length^2 *
                model.equilibrium_distribution(end_to_end_length),
            ZERO,
            ONE,
            POINTS,
        )
        @test normalization - 1.0 <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        normalization = integrate(
            end_to_end_length ->
                model.equilibrium_radial_distribution(end_to_end_length),
            ZERO,
            ONE,
            POINTS,
        )
        @test normalization - 1.0 <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::nondimensional_equilibrium_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        normalization = integrate(
            nondimensional_end_to_end_length_per_link ->
                4.0 *
                pi *
                nondimensional_end_to_end_length_per_link^2 *
                model.nondimensional_equilibrium_distribution(
                    nondimensional_end_to_end_length_per_link,
                ),
            ZERO,
            ONE,
            POINTS,
        )
        @test normalization - 1.0 <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::nondimensional_equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        normalization = integrate(
            nondimensional_end_to_end_length_per_link ->
                model.nondimensional_equilibrium_radial_distribution(
                    nondimensional_end_to_end_length_per_link,
                ),
            ZERO,
            ONE,
            POINTS,
        )
        @test normalization - 1.0 <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::per_link::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::per_link::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::per_link::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::per_link::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::relative::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::relative::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::relative::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::relative::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::relative_helmholtz_free_energy" begin
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
        relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(
            ZERO * number_of_links * link_length,
            temperature,
        )
        @test abs(relative_helmholtz_free_energy_0) <=
              ZERO * number_of_links * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::relative_helmholtz_free_energy_per_link" begin
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
        relative_helmholtz_free_energy_per_link_0 =
            model.relative_helmholtz_free_energy_per_link(
                ZERO * number_of_links * link_length,
                temperature,
            )
        @test abs(relative_helmholtz_free_energy_per_link_0) <=
              ZERO * BOLTZMANN_CONSTANT * temperature
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_relative_helmholtz_free_energy_0 =
            model.nondimensional_relative_helmholtz_free_energy(ZERO)
        @test abs(nondimensional_relative_helmholtz_free_energy_0) <= ZERO * number_of_links
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_relative_helmholtz_free_energy_per_link_0 =
            model.nondimensional_relative_helmholtz_free_energy_per_link(ZERO)
        @test abs(nondimensional_relative_helmholtz_free_energy_per_link_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        equilibrium_radial_distribution_0 =
            model.equilibrium_radial_distribution(ZERO * number_of_links * link_length)
        @test abs(equilibrium_radial_distribution_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::zero::nondimensional_equilibrium_radial_distribution" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_equilibrium_radial_distribution_0 =
            model.equilibrium_radial_distribution(ZERO)
        @test abs(nondimensional_equilibrium_radial_distribution_0) <= ZERO
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::connection::force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::connection::nondimensional_force" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::connection::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::connection::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
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

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        force = model.force(end_to_end_length, temperature)
        helmholtz_free_energy = model.helmholtz_free_energy(end_to_end_length, temperature)
        helmholtz_free_energy_legendre =
            model.legendre.gibbs_free_energy(end_to_end_length, temperature) +
            force * end_to_end_length
        residual_abs = helmholtz_free_energy - helmholtz_free_energy_legendre
        residual_rel = residual_abs / helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        end_to_end_length_per_link =
            nondimensional_end_to_end_length_per_link * link_length
        force = model.force(end_to_end_length, temperature)
        helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        helmholtz_free_energy_per_link_legendre =
            model.legendre.gibbs_free_energy_per_link(end_to_end_length, temperature) +
            force * end_to_end_length_per_link
        residual_abs = helmholtz_free_energy_per_link - helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        force = model.force(end_to_end_length, temperature)
        force_0 = model.force(ZERO*number_of_links*link_length, temperature)
        relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        relative_helmholtz_free_energy_legendre =
            model.legendre.relative_gibbs_free_energy(end_to_end_length, temperature) +
            force * end_to_end_length - force_0*ZERO*number_of_links*link_length
        residual_abs = relative_helmholtz_free_energy - relative_helmholtz_free_energy_legendre
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        end_to_end_length_per_link =
            nondimensional_end_to_end_length_per_link * link_length
        force = model.force(end_to_end_length, temperature)
        force_0 = model.force(ZERO*number_of_links*link_length, temperature)
        relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
        relative_helmholtz_free_energy_per_link_legendre =
            model.legendre.relative_gibbs_free_energy_per_link(end_to_end_length, temperature) +
            force * end_to_end_length_per_link - force_0*ZERO*link_length
        residual_abs = relative_helmholtz_free_energy_per_link - relative_helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::nondimensional_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        nondimensional_force = model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(nondimensional_end_to_end_length_per_link, temperature)
        nondimensional_helmholtz_free_energy_legendre =
            model.legendre.nondimensional_gibbs_free_energy(nondimensional_end_to_end_length_per_link, temperature) +
            nondimensional_force * nondimensional_end_to_end_length
        residual_abs = nondimensional_helmholtz_free_energy - nondimensional_helmholtz_free_energy_legendre
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::nondimensional_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        nondimensional_force = model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link, temperature)
        nondimensional_helmholtz_free_energy_per_link_legendre =
            model.legendre.nondimensional_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link, temperature) +
            nondimensional_force * nondimensional_end_to_end_length_per_link
        residual_abs = nondimensional_helmholtz_free_energy_per_link- nondimensional_helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::nondimensional_relative_helmholtz_free_energy" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        nondimensional_force = model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_force_0 = model.nondimensional_force(ZERO)
        nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(nondimensional_end_to_end_length_per_link)
        nondimensional_relative_helmholtz_free_energy_legendre =
            model.legendre.nondimensional_relative_gibbs_free_energy(nondimensional_end_to_end_length_per_link) +
            nondimensional_force * nondimensional_end_to_end_length - nondimensional_force_0 * ZERO * number_of_links
        residual_abs = nondimensional_relative_helmholtz_free_energy - nondimensional_relative_helmholtz_free_energy_legendre
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre::nondimensional_relative_helmholtz_free_energy_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        nondimensional_force = model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_force_0 = model.nondimensional_force(ZERO)
        nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)
        nondimensional_relative_helmholtz_free_energy_per_link_legendre =
            model.legendre.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_end_to_end_length_per_link) +
            nondimensional_force * nondimensional_end_to_end_length_per_link - nondimensional_force_0 * ZERO
        residual_abs = nondimensional_relative_helmholtz_free_energy_per_link - nondimensional_relative_helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_abs) <= parameters.abs_tol &&
              abs(residual_rel) <= parameters.rel_tol
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre_connection::end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        h = parameters.rel_tol * number_of_links * link_length
        end_to_end_length_from_derivative =
            -(
                model.legendre.relative_gibbs_free_energy(
                    end_to_end_length + 0.5 * h,
                    temperature,
                ) - model.legendre.relative_gibbs_free_energy(
                    end_to_end_length - 0.5 * h,
                    temperature,
                )
            ) / (
                model.force(end_to_end_length + 0.5 * h, temperature) -
                model.force(end_to_end_length - 0.5 * h, temperature)
            )
        residual_abs = end_to_end_length - end_to_end_length_from_derivative
        residual_rel = residual_abs / end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre_connection::end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        h = parameters.rel_tol * number_of_links * link_length
        end_to_end_length_per_link_from_derivative =
            -(
                model.legendre.relative_gibbs_free_energy_per_link(
                    end_to_end_length + 0.5 * h,
                    temperature,
                ) - model.legendre.relative_gibbs_free_energy_per_link(
                    end_to_end_length - 0.5 * h,
                    temperature,
                )
            ) / (
                model.force(end_to_end_length + 0.5 * h, temperature) -
                model.force(end_to_end_length - 0.5 * h, temperature)
            )
        residual_abs =
            end_to_end_length_per_link - end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre_connection::nondimensional_end_to_end_length" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
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
        h = parameters.rel_tol
        nondimensional_end_to_end_length_from_derivative =
            -(
                model.legendre.nondimensional_relative_gibbs_free_energy(
                    nondimensional_end_to_end_length_per_link + 0.5 * h,
                ) - model.legendre.nondimensional_relative_gibbs_free_energy(
                    nondimensional_end_to_end_length_per_link - 0.5 * h,
                )
            ) / (
                model.nondimensional_force(
                    nondimensional_end_to_end_length_per_link + 0.5 * h,
                ) - model.nondimensional_force(
                    nondimensional_end_to_end_length_per_link - 0.5 * h,
                )
            )
        residual_abs =
            nondimensional_end_to_end_length -
            nondimensional_end_to_end_length_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::legendre_connection::nondimensional_end_to_end_length_per_link" begin
    for _ = 1:parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference + parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference + parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference +
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        h = parameters.rel_tol
        nondimensional_end_to_end_length_per_link_from_derivative =
            -(
                model.legendre.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link + 0.5 * h,
                ) - model.legendre.nondimensional_relative_gibbs_free_energy_per_link(
                    nondimensional_end_to_end_length_per_link - 0.5 * h,
                )
            ) / (
                model.nondimensional_force(
                    nondimensional_end_to_end_length_per_link + 0.5 * h,
                ) - model.nondimensional_force(
                    nondimensional_end_to_end_length_per_link - 0.5 * h,
                )
            )
        residual_abs =
            nondimensional_end_to_end_length_per_link -
            nondimensional_end_to_end_length_per_link_from_derivative
        residual_rel = residual_abs / nondimensional_end_to_end_length_per_link
        @test abs(residual_rel) <= h
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::force" begin
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
        force = model.force(end_to_end_length, temperature)
        force_legendre = model.legendre.force(end_to_end_length, temperature)
        residual_abs = force - force_legendre
        residual_rel = residual_abs / force
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_force" begin
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
            model.nondimensional_force(nondimensional_end_to_end_length_per_link)
        nondimensional_force_legendre =
            model.legendre.nondimensional_force(nondimensional_end_to_end_length_per_link)
        residual_abs = nondimensional_force - nondimensional_force_legendre
        residual_rel = residual_abs / nondimensional_force
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::helmholtz_free_energy" begin
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
        helmholtz_free_energy = model.helmholtz_free_energy(end_to_end_length, temperature)
        helmholtz_free_energy_legendre =
            model.legendre.helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs = helmholtz_free_energy - helmholtz_free_energy_legendre
        residual_rel = residual_abs / helmholtz_free_energy
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::helmholtz_free_energy_per_link" begin
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
        helmholtz_free_energy_per_link =
            model.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        force_legendre =
            model.legendre.helmholtz_free_energy_per_link(end_to_end_length, temperature)
        residual_abs = helmholtz_free_energy_per_link - force_legendre
        residual_rel = residual_abs / helmholtz_free_energy_per_link
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::relative_helmholtz_free_energy" begin
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
        relative_helmholtz_free_energy =
            model.relative_helmholtz_free_energy(end_to_end_length, temperature)
        relative_helmholtz_free_energy_legendre =
            model.legendre.relative_helmholtz_free_energy(end_to_end_length, temperature)
        residual_abs =
            relative_helmholtz_free_energy - relative_helmholtz_free_energy_legendre
        residual_rel = residual_abs / relative_helmholtz_free_energy
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::relative_helmholtz_free_energy_per_link" begin
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
        relative_helmholtz_free_energy_per_link =
            model.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)
        force_legendre = model.legendre.relative_helmholtz_free_energy_per_link(
            end_to_end_length,
            temperature,
        )
        residual_abs = relative_helmholtz_free_energy_per_link - force_legendre
        residual_rel = residual_abs / relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_helmholtz_free_energy" begin
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
        nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(
            nondimensional_end_to_end_length_per_link,
            temperature,
        )
        nondimensional_helmholtz_free_energy_legendre =
            model.legendre.nondimensional_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy -
            nondimensional_helmholtz_free_energy_legendre
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_helmholtz_free_energy_per_link" begin
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
        nondimensional_helmholtz_free_energy_per_link =
            model.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        nondimensional_helmholtz_free_energy_per_link_legendre =
            model.legendre.nondimensional_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        residual_abs =
            nondimensional_helmholtz_free_energy_per_link -
            nondimensional_helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy" begin
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
        nondimensional_relative_helmholtz_free_energy =
            model.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_legendre =
            model.legendre.nondimensional_relative_helmholtz_free_energy(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy -
            nondimensional_relative_helmholtz_free_energy_legendre
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_relative_helmholtz_free_energy_per_link" begin
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
        nondimensional_relative_helmholtz_free_energy_per_link =
            model.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_relative_helmholtz_free_energy_per_link_legendre =
            model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_relative_helmholtz_free_energy_per_link -
            nondimensional_relative_helmholtz_free_energy_per_link_legendre
        residual_rel = residual_abs / nondimensional_relative_helmholtz_free_energy_per_link
        @test abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::equilibrium_distribution" begin
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
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        equilibrium_distribution = model.equilibrium_distribution(end_to_end_length)
        equilibrium_distribution_legendre =
            model.legendre.equilibrium_distribution(end_to_end_length)
        residual_abs = equilibrium_distribution - equilibrium_distribution_legendre
        residual_rel = residual_abs / equilibrium_distribution
        @test abs(residual_abs) <= 1.0 / sqrt(number_of_links) ||
              abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_equilibrium_distribution" begin
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
        nondimensional_equilibrium_distribution =
            model.nondimensional_equilibrium_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_equilibrium_distribution_legendre =
            model.legendre.nondimensional_equilibrium_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_equilibrium_distribution -
            nondimensional_equilibrium_distribution_legendre
        residual_rel = residual_abs / nondimensional_equilibrium_distribution
        @test abs(residual_abs) <= 1.0 / sqrt(number_of_links) ||
              abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::equilibrium_radial_distribution" begin
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
        end_to_end_length =
            nondimensional_end_to_end_length_per_link * number_of_links * link_length
        equilibrium_radial_distribution =
            model.equilibrium_radial_distribution(end_to_end_length)
        equilibrium_radial_distribution_legendre =
            model.legendre.equilibrium_radial_distribution(end_to_end_length)
        residual_abs =
            equilibrium_radial_distribution - equilibrium_radial_distribution_legendre
        residual_rel = residual_abs / equilibrium_radial_distribution
        @test abs(residual_abs) <= 1.0 / sqrt(number_of_links) ||
              abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::thermodynamic_limit::nondimensional_equilibrium_radial_distribution" begin
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
        nondimensional_equilibrium_radial_distribution =
            model.nondimensional_equilibrium_radial_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        nondimensional_equilibrium_radial_distribution_legendre =
            model.legendre.nondimensional_equilibrium_radial_distribution(
                nondimensional_end_to_end_length_per_link,
            )
        residual_abs =
            nondimensional_equilibrium_radial_distribution -
            nondimensional_equilibrium_radial_distribution_legendre
        residual_rel = residual_abs / nondimensional_equilibrium_radial_distribution
        @test abs(residual_abs) <= 1.0 / sqrt(number_of_links) ||
              abs(residual_rel) <= 1.0 / sqrt(number_of_links)
    end
end
end
