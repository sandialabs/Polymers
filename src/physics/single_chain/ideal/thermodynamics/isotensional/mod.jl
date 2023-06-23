"""
The ideal chain model thermodynamics in the isotensional ensemble.
"""
module Isotensional

using DocStringExtensions
using Polymers_jll
using ......Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the ideal chain model in the isotensional ensemble.

$(FIELDS)
"""
struct IDEAL
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
    The expected end-to-end length ``\\xi`` as a function of the applied force ``f`` and temperature ``T``.
    """
    end_to_end_length::Function
    """
    The expected end-to-end length per link ``\\xi/N_b=\\ell_b\\gamma`` as a function of the applied force ``f`` and temperature ``T``.
    """
    end_to_end_length_per_link::Function
    """
    The expected nondimensional end-to-end length ``N_b\\gamma=\\xi/\\ell_b`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_end_to_end_length::Function
    """
    The expected nondimensional end-to-end length per link ``\\gamma\\equiv\\xi/N_b\\ell_b`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_end_to_end_length_per_link::Function
    """
    The Gibbs free energy ``\\varphi`` as a function of the applied force ``f`` and temperature ``T``.
    """
    gibbs_free_energy::Function
    """
    The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``.
    """
    gibbs_free_energy_per_link::Function
    """
    The relative Gibbs free energy ``\\Delta\\varphi\\equiv\\varphi(f,T)-\\varphi(0,T)`` as a function of the applied force ``f`` and temperature ``T``.
    """
    relative_gibbs_free_energy::Function
    """
    The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``.
    """
    relative_gibbs_free_energy_per_link::Function
    """
    The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``.
    """
    nondimensional_gibbs_free_energy::Function
    """
    The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``.
    """
    nondimensional_gibbs_free_energy_per_link::Function
    """
    The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_relative_gibbs_free_energy::Function
    """
    The nondimensional relative Gibbs free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta``.
    """
    nondimensional_relative_gibbs_free_energy_per_link::Function
end

"""
The expected end-to-end length ``\\xi`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
\\xi(f, T) = -\\frac{\\partial\\varphi}{\\partial f} = \\frac{N_b\\ell_b^2f}{3kT}.
```

$(TYPEDSIGNATURES)
"""
function end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_end_to_end_length,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        force,
        temperature,
    )
end

"""
The expected end-to-end length per link ``\\xi/N_b=\\ell_b\\gamma`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function end_to_end_length_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_end_to_end_length_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64),
            link_length_i,
            force_i,
            temperature_i,
        ),
        link_length,
        force,
        temperature,
    )
end

"""
The expected nondimensional end-to-end length ``N_b\\gamma\\equiv\\xi/\\ell_b`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_force_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_end_to_end_length,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64),
            number_of_links_i,
            nondimensional_force_i,
        ),
        number_of_links,
        nondimensional_force,
    )
end

"""
The expected nondimensional end-to-end length per link ``\\gamma\\equiv\\xi/N_b\\ell_b`` as a function of the applied nondimensional force ``\\eta``,

```math
\\gamma(\\eta) = -\\frac{\\partial\\varrho}{\\partial\\eta} = \\frac{\\eta}{3}.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length_per_link(
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_force_i -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_end_to_end_length_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64,),
            nondimensional_force_i,
        ),
        nondimensional_force,
    )
end

"""
The Gibbs free energy ``\\varphi`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``,

```math
\\varphi(f, T) = -kT\\ln Z(f, T).
```

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, hinge_mass_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        force,
        temperature,
    )
end

"""
The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, hinge_mass_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        force,
        temperature,
    )
end

"""
The relative Gibbs free energy ``\\Delta\\varphi\\equiv\\varphi(f,T)-\\varphi(0,T)`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
\\Delta\\varphi(f, T) = -\\frac{N_b\\ell_b^2f^2}{6kT}.
```

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_relative_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        force,
        temperature,
    )
end

"""
The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_relative_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64),
            link_length_i,
            force_i,
            temperature_i,
        ),
        link_length,
        force,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, hinge_mass_i, nondimensional_force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            nondimensional_force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional relative Gibbs free energy ``\\beta\\Delta\\varphi=N_b\\Delta\\varrho`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_force_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64),
            number_of_links_i,
            nondimensional_force_i,
        ),
        number_of_links,
        nondimensional_force,
    )
end

"""
The nondimensional relative Gibbs free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta``,

```math
\\Delta\\varrho(\\eta) = -\\frac{\\eta^2}{6}.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy_per_link(
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_force_i -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (Float64,),
            nondimensional_force_i,
        ),
        nondimensional_force,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the ideal chain model in the isotensional ensemble.

$(TYPEDSIGNATURES)
"""
function IDEAL(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return IDEAL(
        number_of_links,
        link_length,
        hinge_mass,
        (force, temperature) ->
            end_to_end_length(number_of_links, link_length, force, temperature),
        (force, temperature) -> end_to_end_length_per_link(link_length, force, temperature),
        nondimensional_force ->
            nondimensional_end_to_end_length(number_of_links, nondimensional_force),
        nondimensional_force ->
            nondimensional_end_to_end_length_per_link(nondimensional_force),
        (force, temperature) ->
            gibbs_free_energy(number_of_links, link_length, hinge_mass, force, temperature),
        (force, temperature) ->
            gibbs_free_energy_per_link(link_length, hinge_mass, force, temperature),
        (force, temperature) ->
            relative_gibbs_free_energy(number_of_links, link_length, force, temperature),
        (force, temperature) ->
            relative_gibbs_free_energy_per_link(link_length, force, temperature),
        (nondimensional_force, temperature) -> nondimensional_gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_force,
            temperature,
        ),
        (nondimensional_force, temperature) -> nondimensional_gibbs_free_energy_per_link(
            link_length,
            hinge_mass,
            nondimensional_force,
            temperature,
        ),
        nondimensional_force -> nondimensional_relative_gibbs_free_energy(
            number_of_links,
            nondimensional_force,
        ),
        nondimensional_force ->
            nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force),
    )
end

end
