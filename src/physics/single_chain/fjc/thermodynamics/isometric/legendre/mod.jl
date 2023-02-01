module Legendre

using .......Polymers: PROJECT_ROOT

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    force::Any
    nondimensional_force::Any
    helmholtz_free_energy::Any
    helmholtz_free_energy_per_link::Any
    relative_helmholtz_free_energy::Any
    relative_helmholtz_free_energy_per_link::Any
    nondimensional_helmholtz_free_energy::Any
    nondimensional_helmholtz_free_energy_per_link::Any
    nondimensional_relative_helmholtz_free_energy::Any
    nondimensional_relative_helmholtz_free_energy_per_link::Any
    equilibrium_distribution::Any
    nondimensional_equilibrium_distribution::Any
    equilibrium_radial_distribution::Any
    nondimensional_equilibrium_radial_distribution::Any
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
        fjc.force =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_force,
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
        fjc.nondimensional_force =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_force,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_helmholtz_free_energy,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_helmholtz_free_energy_per_link,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy_per_link,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy_per_link,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_equilibrium_distribution,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_distribution,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_equilibrium_radial_distribution,
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
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_radial_distribution,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.gibbs_free_energy =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_gibbs_free_energy,
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
        fjc.gibbs_free_energy_per_link =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_gibbs_free_energy_per_link,
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
        fjc.relative_gibbs_free_energy =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_gibbs_free_energy,
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
        fjc.relative_gibbs_free_energy_per_link =
            (end_to_end_length, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_gibbs_free_energy_per_link,
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
        fjc.nondimensional_gibbs_free_energy =
            (nondimensional_end_to_end_length_per_link, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy,
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
        fjc.nondimensional_gibbs_free_energy_per_link =
            (nondimensional_end_to_end_length_per_link, temperature) -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy_per_link,
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
        fjc.nondimensional_relative_gibbs_free_energy =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64),
                number_of_links,
                hinge_mass,
                link_length,
                nondimensional_end_to_end_length_per_link,
            )
        fjc.nondimensional_relative_gibbs_free_energy_per_link =
            nondimensional_end_to_end_length_per_link -> ccall(
                (
                    :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy_per_link,
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
