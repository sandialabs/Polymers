using Test
using Polymers.Physics: BOLTZMANN_CONSTANT, PLANCK_CONSTANT
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
    for _ in parameters.number_of_loops
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
    for _ in parameters.number_of_loops
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        @test FJC(
            parameters.number_of_links_minimum,
            link_length,
            parameters.hinge_mass_reference,
        ).link_length == link_length
    end
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::hinge_mass" begin
    for _ in parameters.number_of_loops
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
        @test FJC(
            parameters.number_of_links_minimum,
            parameters.link_length_reference,
            hinge_mass,
        ).hinge_mass == hinge_mass
    end
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::all_parameters" begin
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
        @test all(
            FJC(number_of_links, link_length, hinge_mass).number_of_links ==
            number_of_links &&
            FJC(number_of_links, link_length, hinge_mass).link_length == link_length &&
            FJC(number_of_links, link_length, hinge_mass).hinge_mass == hinge_mass,
        )
    end
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::equilibrium_distribution" begin
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
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
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
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
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
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
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
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
    for _ in parameters.number_of_loops
        number_of_links =
            rand(parameters.number_of_links_minimum:parameters.number_of_links_maximum)
        link_length =
            parameters.link_length_reference .+
            parameters.link_length_scale * (0.5 - rand())
        hinge_mass =
            parameters.hinge_mass_reference .+ parameters.hinge_mass_scale * (0.5 - rand())
        model = FJC(number_of_links, link_length, hinge_mass)
        nondimensional_end_to_end_length_per_link =
            parameters.nondimensional_end_to_end_length_per_link_reference .+
            parameters.nondimensional_end_to_end_length_per_link_scale * (0.5 - rand())
        temperature =
            parameters.temperature_reference .+
            parameters.temperature_scale * (0.5 - rand())
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
