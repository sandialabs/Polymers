module Isotensional

using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    legendre::Any
    end_to_end_length::Any
    end_to_end_length_per_link::Any
    nondimensional_end_to_end_length::Any
    nondimensional_end_to_end_length_per_link::Any
    gibbs_free_energy::Any
    gibbs_free_energy_per_link::Any
    relative_gibbs_free_energy::Any
    relative_gibbs_free_energy_per_link::Any
    nondimensional_gibbs_free_energy::Any
    nondimensional_gibbs_free_energy_per_link::Any
    nondimensional_relative_gibbs_free_energy::Any
    nondimensional_relative_gibbs_free_energy_per_link::Any
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.legendre = Legendre.FJC(number_of_links, link_length, hinge_mass)
        fjc.end_to_end_length =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.end_to_end_length_per_link =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.nondimensional_end_to_end_length =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
            )
        fjc.nondimensional_end_to_end_length_per_link =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
            )
        fjc.gibbs_free_energy =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.gibbs_free_energy_per_link =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_gibbs_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.relative_gibbs_free_energy =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_relative_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.relative_gibbs_free_energy_per_link =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_relative_gibbs_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                force,
                temperature,
            )
        fjc.nondimensional_gibbs_free_energy =
            (nondimensional_force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
                temperature,
            )
        fjc.nondimensional_gibbs_free_energy_per_link =
            (nondimensional_force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_gibbs_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
                temperature,
            )
        fjc.nondimensional_relative_gibbs_free_energy =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
            )
        fjc.nondimensional_relative_gibbs_free_energy_per_link =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
            )
        return fjc
    end
end

end
