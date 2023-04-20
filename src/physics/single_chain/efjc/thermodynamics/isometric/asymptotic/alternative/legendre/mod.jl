"""
The extensible freely-jointed chain (EFJC) model thermodynamics in the isometric ensemble approximated using an alternative asymptotic approach and a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using .........Polymers: PROJECT_ROOT

import ........Physics: BOLTZMANN_CONSTANT

"""
The structure of the thermodynamics of the EFJC model in the isometric ensemble approximated using an alternative asymptotic approach and a Legendre transformation.
$(FIELDS)
"""
struct EFJC
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
    The expected force ``f`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    force::Function
    """
    The expected nondimensional force ``\\eta`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_force::Function
    """
    The Helmholtz free energy ``\\psi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    helmholtz_free_energy::Function
    """
    The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    helmholtz_free_energy_per_link::Function
    """
    The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    relative_helmholtz_free_energy::Function
    """
    The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    relative_helmholtz_free_energy_per_link::Function
    """
    The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy::Function
    """
    The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``.
    """
    nondimensional_helmholtz_free_energy_per_link::Function
    """
    The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_relative_helmholtz_free_energy::Function
    """
    The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``
    """
    nondimensional_relative_helmholtz_free_energy_per_link::Function
end

"""
The expected force as a function ``f`` of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        link_stiffness,
        end_to_end_length,
        temperature,
    )
end

"""
The expected nondimensional force as a function ``\\eta`` of the applied nondimensional end-to-end length per link ``\\gamma``,
parameterized by the link length ``\\ell_b`` and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (nondimensional_link_stiffness_i, nondimensional_end_to_end_length_per_link_i) ->
            ccall(
                (
                    :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_nondimensional_force,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (Float64, Float64),
                nondimensional_link_stiffness_i,
                nondimensional_end_to_end_length_per_link_i,
            ),
        nondimensional_link_stiffness,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The Helmholtz free energy ``\\psi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, link stiffness ``k_0``, and hinge mass ``m``,

```math
\\psi(\\xi, T) \\sim \\varphi\\left[f(\\xi, T)\\right] + \\xi f(\\xi, T) \\quad \\text{for } N_b\\gg 1,
```

where ``f(\\xi, T)`` is given by the Legendre transformation approximation above.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        end_to_end_length,
        temperature,
    )
end

"""
The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, link stiffness ``k_0``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_relative_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        link_stiffness,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and link stiffness ``k_0``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    link_stiffness::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_relative_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            link_stiffness_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        link_stiffness,
        end_to_end_length,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_nondimensional_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_nondimensional_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_link_stiffness,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b`` and nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_nondimensional_relative_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_link_stiffness_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_link_stiffness,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the nondimensional link stiffness ``\\kappa\\equiv\\beta k_0\\ell_b^2``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    nondimensional_link_stiffness::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (nondimensional_link_stiffness_i, nondimensional_end_to_end_length_per_link_i) ->
            ccall(
                (
                    :physics_single_chain_efjc_thermodynamics_isometric_asymptotic_alternative_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
                    string(PROJECT_ROOT, "target/debug/libpolymers"),
                ),
                Float64,
                (Float64, Float64),
                nondimensional_link_stiffness_i,
                nondimensional_end_to_end_length_per_link_i,
            ),
        nondimensional_link_stiffness,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the EFJC model in the isometric ensemble approximated using an alternative asymptotic approach and a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function EFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
)
    BOLTZMANN_CONSTANT::Float64 = 8.314462618
    return EFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        (end_to_end_length, temperature) -> force(
            number_of_links,
            link_length,
            link_stiffness,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) -> nondimensional_force(
            link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
            nondimensional_end_to_end_length_per_link,
        ),
        (end_to_end_length, temperature) -> helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy(
            number_of_links,
            link_length,
            link_stiffness,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            link_stiffness,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_relative_helmholtz_free_energy(
                number_of_links,
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                nondimensional_end_to_end_length_per_link,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_relative_helmholtz_free_energy_per_link(
                link_stiffness * link_length^2 / BOLTZMANN_CONSTANT / temperature,
                nondimensional_end_to_end_length_per_link,
            ),
    )
end

end
