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
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    @test all(
        map(
            number_of_links_i ->
                FJC(
                    number_of_links_i,
                    parameters.link_length_reference,
                    parameters.hinge_mass_reference,
                ).number_of_links,
            number_of_links_list,
        ) == number_of_links_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::link_length" begin
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    @test all(
        map(
            link_length_i ->
                FJC(
                    parameters.number_of_links_minimum,
                    link_length_i,
                    parameters.hinge_mass_reference,
                ).link_length,
            link_length_list,
        ) == link_length_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::hinge_mass" begin
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    @test all(
        map(
            hinge_mass_i ->
                FJC(
                    parameters.number_of_links_minimum,
                    parameters.link_length_reference,
                    hinge_mass_i,
                ).hinge_mass,
            hinge_mass_list,
        ) == hinge_mass_list,
    )
end
@testset "physics::single_chain::fjc::thermodynamics::isometric::test::base::all_parameters" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                FJC(number_of_links_i, link_length_i, hinge_mass_i).number_of_links,
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) == number_of_links_list .&&
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                FJC(number_of_links_i, link_length_i, hinge_mass_i).link_length,
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) == link_length_list .&&
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                FJC(number_of_links_i, link_length_i, hinge_mass_i).hinge_mass,
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) == hinge_mass_list,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::equilibrium_distribution" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    normalization =
        (number_of_links_i, link_length_i, hinge_mass_i) -> integrate(
            end_to_end_length ->
                4.0 *
                pi *
                end_to_end_length^2 *
                FJC(number_of_links_i, link_length_i, hinge_mass_i).equilibrium_distribution(
                    end_to_end_length,
                ),
            ZERO * number_of_links_i * link_length_i,
            ONE * number_of_links_i * link_length_i,
            POINTS,
        )
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                abs(normalization(number_of_links_i, link_length_i, hinge_mass_i) - 1.0),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) .<= parameters.rel_tol,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::equilibrium_radial_distribution" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    normalization =
        (number_of_links_i, link_length_i, hinge_mass_i) -> integrate(
            end_to_end_length -> FJC(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
            ).equilibrium_radial_distribution(
                end_to_end_length,
            ),
            ZERO * number_of_links_i * link_length_i,
            ONE * number_of_links_i * link_length_i,
            POINTS,
        )
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                abs(normalization(number_of_links_i, link_length_i, hinge_mass_i) - 1.0),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) .<= parameters.rel_tol,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::nondimensional_equilibrium_distribution" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    normalization =
        (number_of_links_i, link_length_i, hinge_mass_i) -> integrate(
            nondimensional_end_to_end_length_per_link ->
                4.0 *
                pi *
                nondimensional_end_to_end_length_per_link^2 *
                FJC(
                    number_of_links_i,
                    link_length_i,
                    hinge_mass_i,
                ).nondimensional_equilibrium_distribution(
                    nondimensional_end_to_end_length_per_link,
                ),
            ZERO,
            ONE,
            POINTS,
        )
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                abs(normalization(number_of_links_i, link_length_i, hinge_mass_i) - 1.0),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) .<= parameters.rel_tol,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::normalization::nondimensional_equilibrium_radial_distribution" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    normalization =
        (number_of_links_i, link_length_i, hinge_mass_i) -> integrate(
            nondimensional_end_to_end_length_per_link -> FJC(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
            ).nondimensional_equilibrium_radial_distribution(
                nondimensional_end_to_end_length_per_link,
            ),
            ZERO,
            ONE,
            POINTS,
        )
    @test all(
        map(
            (number_of_links_i, link_length_i, hinge_mass_i) ->
                abs(normalization(number_of_links_i, link_length_i, hinge_mass_i) - 1.0),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
        ) .<= parameters.rel_tol,
    )
end

@testset "physics::single_chain::fjc::thermodynamics::isometric::test::nondimensional::force" begin
    number_of_links_list = rand(
        parameters.number_of_links_minimum:parameters.number_of_links_maximum,
        parameters.number_of_loops,
    )
    link_length_list =
        parameters.link_length_reference .+
        parameters.link_length_scale * (0.5 .- rand(parameters.number_of_loops))
    hinge_mass_list =
        parameters.hinge_mass_reference .+
        parameters.hinge_mass_scale * (0.5 .- rand(parameters.number_of_loops))
    temperature_list =
        parameters.temperature_reference .+
        parameters.temperature_scale * (0.5 .- rand(parameters.number_of_loops))
    nondimensional_end_to_end_length_per_link_list =
        parameters.nondimensional_end_to_end_length_per_link_reference .+
        parameters.nondimensional_end_to_end_length_per_link_scale *
        (0.5 .- rand(parameters.number_of_loops))
    nondimensional_force =
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> FJC(number_of_links_i, link_length_i, hinge_mass_i).nondimensional_force(
            nondimensional_end_to_end_length_per_link_i,
        )
    end_to_end_length_list =
        nondimensional_end_to_end_length_per_link_list .* number_of_links_list .*
        link_length_list
    force =
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ) -> FJC(number_of_links_i, link_length_i, hinge_mass_i).force(
            end_to_end_length_i,
            temperature_i,
        )
    residual_abs =
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            end_to_end_length_i,
            temperature_i,
        ) ->
            force(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                end_to_end_length_i,
                temperature_i,
            ) / BOLTZMANN_CONSTANT / temperature_i * link_length_i - nondimensional_force(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                nondimensional_end_to_end_length_per_link_i,
            )
    residual_rel =
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            end_to_end_length_i,
            temperature_i,
        ) ->
            residual_abs(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                nondimensional_end_to_end_length_per_link_i,
                end_to_end_length_i,
                temperature_i,
            ) / nondimensional_force(
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                nondimensional_end_to_end_length_per_link_i,
            )
    @test all(
        map(
            (
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                nondimensional_end_to_end_length_per_link_i,
                end_to_end_length_i,
                temperature_i,
            ) -> abs(
                residual_abs(
                    number_of_links_i,
                    link_length_i,
                    hinge_mass_i,
                    nondimensional_end_to_end_length_per_link_i,
                    end_to_end_length_i,
                    temperature_i,
                ),
            ),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
            nondimensional_end_to_end_length_per_link_list,
            end_to_end_length_list,
            temperature_list,
        ) .<= parameters.abs_tol .&&
        map(
            (
                number_of_links_i,
                link_length_i,
                hinge_mass_i,
                nondimensional_end_to_end_length_per_link_i,
                end_to_end_length_i,
                temperature_i,
            ) -> abs(
                residual_rel(
                    number_of_links_i,
                    link_length_i,
                    hinge_mass_i,
                    nondimensional_end_to_end_length_per_link_i,
                    end_to_end_length_i,
                    temperature_i,
                ),
            ),
            number_of_links_list,
            link_length_list,
            hinge_mass_list,
            nondimensional_end_to_end_length_per_link_list,
            end_to_end_length_list,
            temperature_list,
        ) .<= parameters.rel_tol,
    )
end
