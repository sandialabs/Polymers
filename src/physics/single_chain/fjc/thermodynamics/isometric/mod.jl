"""
Module docstring!
"""
module Isometric

using DocStringExtensions
using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

"""
    test

This is a test.

$(FIELDS)
"""
mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    legendre::Any
    """
    The force ``f`` as a function of the end-to-end length ``\\xi`` and the temperature ``T``.
    """
    force::Function
    nondimensional_force::Function
    helmholtz_free_energy::Function
    helmholtz_free_energy_per_link::Function
    relative_helmholtz_free_energy::Function
    relative_helmholtz_free_energy_per_link::Function
    nondimensional_helmholtz_free_energy::Function
    nondimensional_helmholtz_free_energy_per_link::Function
    nondimensional_relative_helmholtz_free_energy::Function
    nondimensional_relative_helmholtz_free_energy_per_link::Function
    equilibrium_distribution::Function
    nondimensional_equilibrium_distribution::Function
    equilibrium_radial_distribution::Function
    """
        test also

    This is also a test.

    $(SIGNATURES)
    """
    nondimensional_equilibrium_radial_distribution::Any
    """
        test also also

    This is also also a test.
    """
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.legendre = Legendre.FJC(number_of_links, link_length, hinge_mass)
        fjc.force = function (
            end_to_end_length::Union{Float64,Union{Vector,Matrix}},
            temperature::Union{Float64,Union{Vector,Matrix}},
        )
            if isa(end_to_end_length, Float64) && isa(temperature, Float64)
                return ccall(
                    (
                        :physics_single_chain_fjc_thermodynamics_isometric_force,
                        string(PROJECT_ROOT, "target/debug/libpolymers"),
                    ),
                    Float64,
                    (UInt8, Float64, Float64, Float64, Float64),
                    number_of_links,
                    hinge_mass,
                    link_length,
                    end_to_end_length,
                    temperature,
                )
            elseif isa(end_to_end_length, Union{Vector,Matrix}) && isa(temperature, Float64)
                return broadcast(
                    end_to_end_length_i -> ccall(
                        (
                            :physics_single_chain_fjc_thermodynamics_isometric_force,
                            string(PROJECT_ROOT, "target/debug/libpolymers"),
                        ),
                        Float64,
                        (UInt8, Float64, Float64, Float64, Float64),
                        number_of_links,
                        hinge_mass,
                        link_length,
                        end_to_end_length_i,
                        temperature,
                    ),
                    end_to_end_length,
                )
            elseif isa(end_to_end_length, Float64) && isa(temperature, Union{Vector,Matrix})
                return broadcast(
                    temperature_i -> ccall(
                        (
                            :physics_single_chain_fjc_thermodynamics_isometric_force,
                            string(PROJECT_ROOT, "target/debug/libpolymers"),
                        ),
                        Float64,
                        (UInt8, Float64, Float64, Float64, Float64),
                        number_of_links,
                        hinge_mass,
                        link_length,
                        end_to_end_length,
                        temperature_i,
                    ),
                    temperature,
                )
            elseif isa(end_to_end_length, Union{Vector,Matrix}) &&
                   isa(temperature, Union{Vector,Matrix})
                return broadcast(
                    (end_to_end_length_i, temperature_i) -> ccall(
                        (
                            :physics_single_chain_fjc_thermodynamics_isometric_force,
                            string(PROJECT_ROOT, "target/debug/libpolymers"),
                        ),
                        Float64,
                        (UInt8, Float64, Float64, Float64, Float64),
                        number_of_links,
                        hinge_mass,
                        link_length,
                        end_to_end_length_i,
                        temperature_i,
                    ),
                    end_to_end_length,
                    temperature,
                )
            end
        end
        fjc.nondimensional_force =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_force,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.helmholtz_free_energy =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_helmholtz_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
                temperature,
            )
        fjc.helmholtz_free_energy_per_link =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_helmholtz_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
                temperature,
            )
        fjc.relative_helmholtz_free_energy =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_relative_helmholtz_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
                temperature,
            )
        fjc.relative_helmholtz_free_energy_per_link =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_relative_helmholtz_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
                temperature,
            )
        fjc.nondimensional_helmholtz_free_energy =
            (nondimensional_end_to_end_length_per_link, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_helmholtz_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        fjc.nondimensional_helmholtz_free_energy_per_link =
            (nondimensional_end_to_end_length_per_link, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_helmholtz_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
                temperature,
            )
        fjc.nondimensional_relative_helmholtz_free_energy =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.nondimensional_relative_helmholtz_free_energy_per_link =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.equilibrium_distribution =
            end_to_end_length -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_equilibrium_distribution,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
            )
        fjc.nondimensional_equilibrium_distribution =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_equilibrium_distribution,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.equilibrium_radial_distribution =
            end_to_end_length -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_equilibrium_radial_distribution,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                end_to_end_length,
            )
        fjc.nondimensional_equilibrium_radial_distribution =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_equilibrium_radial_distribution,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        return fjc
    end
end

end
