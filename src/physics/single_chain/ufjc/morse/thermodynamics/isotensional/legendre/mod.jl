"""
The Morse potential freely-jointed chain (Morse-FJC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using Polymers_jll
using ........Polymers: PROJECT_ROOT

import .......Physics: BOLTZMANN_CONSTANT

"""
The structure of the thermodynamics of the Morse-FJC model in the isotensional ensemble approximated using a Legendre transformation.
$(FIELDS)
"""
struct MORSEFJC
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
    The stiffness of each link in the chain ``k_0`` in units of J/(mol⋅nm^2).
    """
    link_stiffness::Float64
    """
    The energy of each link in the chain ``u_0`` in units of J/mol.
    """
    link_energy::Float64
    """
    The Helmholtz free energy ``\\psi`` as a function of the applied force ``f`` and temperature ``T``.
    """
    helmholtz_free_energy::Function
    """
    The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied force ``f`` and temperature ``T``.
    """
    helmholtz_free_energy_per_link::Function
    """
    The relative helmholtz free energy ``\\Delta\\psi\\equiv\\psi(f,T)-\\psi(0,T)`` as a function of the applied force ``f`` and temperature ``T``.
    """
    relative_helmholtz_free_energy::Function
    """
    The relative helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied force ``f`` and temperature ``T``.
    """
    relative_helmholtz_free_energy_per_link::Function
    """
    The nondimensional helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy::Function
    """
    The nondimensional helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy_per_link::Function
    """
    The nondimensional relative helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_relative_helmholtz_free_energy::Function
    """
    The nondimensional relative helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_relative_helmholtz_free_energy_per_link::Function
end

"""
The Helmholtz free energy ``\\psi`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and link stiffness ``k_0``.

```math
\\psi(f, T) \\sim \\varphi(f, T) + f \\xi(f, T) \\quad \\text{for } N_b\\gg 1.
```

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    link_energy::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        link_energy,
        force,
        temperature,
    )
end

"""
The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b``, hinge mass ``m``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    link_energy::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        link_stiffness,
        link_energy,
        force,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(f,T)-\\psi(0,T)`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    link_energy::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        link_stiffness,
        link_energy,
        force,
        temperature,
    )
end

"""
The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    link_energy::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, link_stiffness_i, link_energy_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64, Float64, Float64),
            link_length_i,
            link_stiffness_i,
            link_energy_i,
            force_i,
            temperature_i,
        ),
        link_length,
        link_stiffness,
        link_energy,
        force,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_energy::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_link_energy,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the link length ``\\ell_b``, hinge mass ``m``, and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_energy::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_link_energy,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the number of links ``N_b`` and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_energy::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
        ),
        number_of_links,
        nondimensional_link_stiffness,
        nondimensional_link_energy,
        nondimensional_force,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional force ``\\eta``
parameterized by the nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_energy::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64),
            nondimensional_link_stiffness_i,
            nondimensional_link_energy_i,
            nondimensional_force_i,
        ),
        nondimensional_link_stiffness,
        nondimensional_link_energy,
        nondimensional_force,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the Morse-FJC model in the isotensional ensemble approximated using a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function MORSEFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
    link_energy::Float64,
)
    BOLTZMANN_CONSTANT::Float64 = 8.314462618
    return MORSEFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        link_energy,
        (force, temperature) -> helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
            force,
            temperature,
        ),
        (force, temperature) -> helmholtz_free_energy_per_link(
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
            force,
            temperature,
        ),
        (force, temperature) -> relative_helmholtz_free_energy(
            number_of_links,
            link_length,
            link_stiffness,
            link_energy,
            force,
            temperature,
        ),
        (force, temperature) -> relative_helmholtz_free_energy_per_link(
            link_length,
            link_stiffness,
            link_energy,
            force,
            temperature,
        ),
        (nondimensional_force, temperature) -> nondimensional_helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            link_energy / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
            temperature,
        ),
        (nondimensional_force, temperature) ->
            nondimensional_helmholtz_free_energy_per_link(
                link_length,
                hinge_mass,
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                link_energy / BOLTZMANN_CONSTANT / temperature,
                nondimensional_force,
                temperature,
            ),
        (nondimensional_force, temperature) ->
            nondimensional_relative_helmholtz_free_energy(
                number_of_links,
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                link_energy / BOLTZMANN_CONSTANT / temperature,
                nondimensional_force,
            ),
        (nondimensional_force, temperature) ->
            nondimensional_relative_helmholtz_free_energy_per_link(
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                link_energy / BOLTZMANN_CONSTANT / temperature,
                nondimensional_force,
            ),
    )
end

end
