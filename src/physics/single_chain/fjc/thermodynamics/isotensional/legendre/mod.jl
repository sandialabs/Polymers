module Legendre

using .......Polymers: PROJECT_ROOT

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    helmholtz_free_energy::Any
    helmholtz_free_energy_per_link::Any
    relative_helmholtz_free_energy::Any
    relative_helmholtz_free_energy_per_link::Any
    nondimensional_helmholtz_free_energy::Any
    nondimensional_helmholtz_free_energy_per_link::Any
    nondimensional_relative_helmholtz_free_energy::Any
    nondimensional_relative_helmholtz_free_energy_per_link::Any
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.helmholtz_free_energy =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_helmholtz_free_energy,
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
        fjc.helmholtz_free_energy_per_link =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_helmholtz_free_energy_per_link,
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
        fjc.relative_helmholtz_free_energy =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy,
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
        fjc.relative_helmholtz_free_energy_per_link =
            (force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy_per_link,
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
        fjc.nondimensional_helmholtz_free_energy =
            (nondimensional_force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy,
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
        fjc.nondimensional_helmholtz_free_energy_per_link =
            (nondimensional_force, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy_per_link,
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
        fjc.nondimensional_relative_helmholtz_free_energy =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_force,
            )
        fjc.nondimensional_relative_helmholtz_free_energy_per_link =
            nondimensional_force -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
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
