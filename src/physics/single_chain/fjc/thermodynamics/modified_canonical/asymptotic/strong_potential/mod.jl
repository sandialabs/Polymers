"""
The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
"""
module StrongPotential

using DocStringExtensions
using Polymers_jll
using ........Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.

$(FIELDS)
"""
struct FJC
    """
    The number of links in the chain ``N_b``.
    """
    number_of_links::UInt8
    """
    The length of each link in the chain ``\\ell_b`` in units of nm.
    """
    link_length::Float64
    """
    The mass of each hinge in the chain ``m`` in units of kg/mol.
    """
    hinge_mass::Float64
    """
    The expected force ``f`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    force::Function
    """
    The expected nondimensional force ``\\eta`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_force::Function
    """
    The Helmholtz free energy ``\\psi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    helmholtz_free_energy::Function
    """
    The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    helmholtz_free_energy_per_link::Function
    """
    The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    relative_helmholtz_free_energy::Function
    """
    The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    relative_helmholtz_free_energy_per_link::Function
    """
    The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy::Function
    """
    The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy_per_link::Function
    """
    The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_relative_helmholtz_free_energy::Function
    """
    The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_relative_helmholtz_free_energy_per_link::Function
end

"""
The expected force ``f`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_force,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The expected nondimensional force ``\\eta`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``, given by [Buche and Rimsza](https://doi.org/10.48550/arXiv.2309.01009) as

```math
\\eta(\\gamma) = \\eta_0(\\gamma) - \\frac{1}{N_b\\varpi}\\left[\\eta_0(\\gamma)\\eta_0'(\\gamma) - \\frac{\\eta_0''(\\gamma)}{2N_b}\\right],
```

where ``\\eta_0(\\gamma)`` is the isometric mechanical response.

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_nondimensional_force,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The Helmholtz free energy ``\\psi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_relative_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_nondimensional_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_nondimensional_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_nondimensional_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_strong_potential_nondimensional_relative_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        (potential_distance, potential_stiffness, temperature) -> force(
            number_of_links,
            link_length,
            potential_distance,
            potential_stiffness,
            temperature,
        ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_force(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (potential_distance, potential_stiffness, temperature) -> helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            potential_distance,
            potential_stiffness,
            temperature,
        ),
        (potential_distance, potential_stiffness, temperature) ->
            helmholtz_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (potential_distance, potential_stiffness, temperature) ->
            relative_helmholtz_free_energy(
                number_of_links,
                link_length,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (potential_distance, potential_stiffness, temperature) ->
            relative_helmholtz_free_energy_per_link(
                number_of_links,
                link_length,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ) -> nondimensional_helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ),
        (
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ) -> nondimensional_helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_relative_helmholtz_free_energy(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_relative_helmholtz_free_energy_per_link(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
    )
end

end
