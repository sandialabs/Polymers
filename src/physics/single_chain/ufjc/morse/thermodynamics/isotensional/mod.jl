"""
The Morse potential freely-jointed chain (Morse-FJC) model thermodynamics in the isotensional ensemble.
"""
module Isotensional

using DocStringExtensions
using .......Polymers: PROJECT_ROOT

import ......Physics: BOLTZMANN_CONSTANT

include("asymptotic/mod.jl")
include("legendre/mod.jl")

"""
The structure of the thermodynamics of the Morse-FJC model in the isotensional ensemble.

$(FIELDS)
"""
struct MORSE
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
    The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    """
    asymptotic::Any
    """
    The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    """
    legendre::Any
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
The expected end-to-end length ``N_b\\gamma=\\xi/\\ell_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``,

```math
\\xi(f, T) = -\\frac{\\partial\\varphi}{\\partial f}.
```

$(TYPEDSIGNATURES)
"""
function end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, link_stiffness_i, force_i, temperature_i) ->
            ccall(
                (
                    :physics_single_chain_ufjc_morse_thermodynamics_isotensional_end_to_end_length,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links_i,
                link_length_i,
                link_stiffness_i,
                force_i,
                temperature_i,
            ),
        number_of_links,
        link_length,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The expected end-to-end length per link ``\\xi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function end_to_end_length_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, link_stiffness_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_end_to_end_length_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64, Float64, Float64),
            link_length_i,
            link_stiffness_i,
            force_i,
            temperature_i,
        ),
        link_length,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The expected nondimensional end-to-end length ``N_b\\gamma=\\xi/\\ell_b`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the number of links ``N_b`` and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_link_stiffness_i, nondimensional_force_i) ->
            ccall(
                (
                    :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_end_to_end_length,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64),
                number_of_links_i,
                nondimensional_link_stiffness_i,
                nondimensional_force_i,
            ),
        number_of_links,
        nondimensional_link_stiffness,
        nondimensional_force,
    )
end

"""
The expected nondimensional end-to-end length per link ``\\gamma\\equiv \\xi/N_b\\ell_b`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``,
given by [Buche et al.](https://doi.org/10.1103/PhysRevE.106.024502) as

```math
\\gamma(\\eta) = -\\frac{\\partial}{\\partial\\eta}\\,\\ln\\left[\\int \\frac{\\sinh(s\\eta)}{s\\eta}\\,e^{-\\beta u(s)}s^2\\,ds\\right].
```

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length_per_link(
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (nondimensional_link_stiffness_i, nondimensional_force_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_end_to_end_length_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64),
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
        ),
        nondimensional_link_stiffness,
        nondimensional_force,
    )
end

"""
The Gibbs free energy ``\\varphi`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and link stiffness ``k_0``.

```math
\\varphi(f, T) = -kT\\ln Z(f, T).
```

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b``, hinge mass ``m``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, hinge_mass_i, link_stiffness_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The relative Gibbs free energy ``\\Delta\\varphi\\equiv\\varphi(f,T)-\\varphi(0,T)`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, link_stiffness_i, force_i, temperature_i) ->
            ccall(
                (
                    :physics_single_chain_ufjc_morse_thermodynamics_isotensional_relative_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64, Float64, Float64),
                number_of_links_i,
                link_length_i,
                link_stiffness_i,
                force_i,
                temperature_i,
            ),
        number_of_links,
        link_length,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, link_stiffness_i, force_i, temperature_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64, Float64, Float64),
            link_length_i,
            link_stiffness_i,
            force_i,
            temperature_i,
        ),
        link_length,
        link_stiffness,
        force,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta`` and temperature ``T``,
parameterized by the link length ``\\ell_b``, hinge mass ``m``, and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_force,
        temperature,
    )
end

"""
The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the number of links ``N_b`` and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_link_stiffness_i, nondimensional_force_i) ->
            ccall(
                (
                    :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (UInt8, Float64, Float64),
                number_of_links_i,
                nondimensional_link_stiffness_i,
                nondimensional_force_i,
            ),
        number_of_links,
        nondimensional_link_stiffness,
        nondimensional_force,
    )
end

"""
The nondimensional relative Gibbs free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional force ``\\eta``,
parameterized by the nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``,
given by [Buche et al.](https://doi.org/10.1103/PhysRevE.106.024502) as

```math
\\Delta\\varrho(\\eta) = \\ln\\left[\\int \\frac{\\sinh(s\\eta)}{s\\eta}\\,e^{-\\beta u(s)}s^2\\,ds\\right] - \\ln\\left[\\int e^{-\\beta u(s)}s^2\\,ds\\right],
```

where the nondimensional link potential ``\\beta u`` is given by [Morse](https://doi.org/10.1103/PhysRev.34.57) as

```math
\\beta u(\\lambda) = \\varepsilon\\left[1 - e^{\\alpha(\\lambda - 1)}\\right]^2,
```

where ``\\varepsilon\\equiv\\beta u_b=\\kappa/72`` is the nondimensional potential energy scale, ``\\alpha\\equiv a\\ell_b=\\sqrt{\\kappa/2\\varepsilon}`` is the nondimensional Morse parameter, ``\\kappa\\equiv\\beta k_b\\ell_b^2`` is the nondimensional link stiffness, and ``\\lambda\\equiv\\ell/\\ell_b`` is the nondimensional link stretch.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy_per_link(
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_force::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (nondimensional_link_stiffness_i, nondimensional_force_i) -> ccall(
            (
                :physics_single_chain_ufjc_morse_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64),
            nondimensional_link_stiffness_i,
            nondimensional_force_i,
        ),
        nondimensional_link_stiffness,
        nondimensional_force,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the Morse-FJC model in the isotensional ensemble.

$(TYPEDSIGNATURES)
"""
function MORSE(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
    link_energy::Float64,
)
    BOLTZMANN_CONSTANT::Float64 = 8.314462618
    return MORSE(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        link_energy,
        Asymptotic.MORSE(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
        ),
        Legendre.MORSE(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
        ),
        (force, temperature) -> end_to_end_length(
            number_of_links,
            link_length,
            link_stiffness,
            force,
            temperature,
        ),
        (force, temperature) ->
            end_to_end_length_per_link(link_length, link_stiffness, force, temperature),
        (nondimensional_force, temperature) -> nondimensional_end_to_end_length(
            number_of_links,
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
        ),
        (nondimensional_force, temperature) -> nondimensional_end_to_end_length_per_link(
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
        ),
        (force, temperature) -> gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            force,
            temperature,
        ),
        (force, temperature) -> gibbs_free_energy_per_link(
            link_length,
            hinge_mass,
            link_stiffness,
            force,
            temperature,
        ),
        (force, temperature) -> relative_gibbs_free_energy(
            number_of_links,
            link_length,
            link_stiffness,
            force,
            temperature,
        ),
        (force, temperature) -> relative_gibbs_free_energy_per_link(
            link_length,
            link_stiffness,
            force,
            temperature,
        ),
        (nondimensional_force, temperature) -> nondimensional_gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
            temperature,
        ),
        (nondimensional_force, temperature) -> nondimensional_gibbs_free_energy_per_link(
            link_length,
            hinge_mass,
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
            temperature,
        ),
        (nondimensional_force, temperature) -> nondimensional_relative_gibbs_free_energy(
            number_of_links,
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_force,
        ),
        (nondimensional_force, temperature) ->
            nondimensional_relative_gibbs_free_energy_per_link(
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                nondimensional_force,
            ),
    )
end

end
