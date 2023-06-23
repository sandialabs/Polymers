"""
The worm-like chain (WLC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using Polymers_jll
using .......Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the WLC model in the isometric ensemble approximated using a Legendre transformation.
$(FIELDS)
"""
struct WLC
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
    The persistance length of the chain in units of nm.
    """
    persistance_length::Float64
    nondimensional_persistance_length::Float64
    """
    The Gibbs free energy ``\\varphi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    gibbs_free_energy::Function
    """
    The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    gibbs_free_energy_per_link::Function
    """
    The relative Gibbs free energy ``\\Delta\\varphi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    relative_gibbs_free_energy::Function
    """
    The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    relative_gibbs_free_energy_per_link::Function
    """
    The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``.
    """
    nondimensional_gibbs_free_energy::Function
    """
    The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``.
    """
    nondimensional_gibbs_free_energy_per_link::Function
    """
    The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_relative_gibbs_free_energy::Function
    """
    The nondimensional relative Gibbs free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``
    """
    nondimensional_relative_gibbs_free_energy_per_link::Function
end

"""
The Gibbs free energy ``\\varphi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and persistance length ``\\ell_p``,

```math
\\varphi(\\xi, T) \\sim \\psi(\\xi, T) - \\xi f(\\xi, T) \\quad \\text{for } N_b\\gg 1,
```

where ``f(\\xi, T)`` is given by the Legendre transformation approximation above.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    persistance_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        persistance_length,
        end_to_end_length,
        temperature,
    )
end

"""
The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and persistance length ``\\ell_p``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    persistance_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        persistance_length,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\varphi\\equiv\\varphi(\\xi,T)-\\varphi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and persistance length ``\\ell_p``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    persistance_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_relative_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        persistance_length,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and persistance length ``\\ell_p``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    persistance_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_relative_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            persistance_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        persistance_length,
        end_to_end_length,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and nondimensional persistance length ``\\zeta``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_persistance_length::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_persistance_length,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, hinge mass ``m``, and nondimensional persistance length ``\\zeta``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_persistance_length::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_persistance_length,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b`` and nondimensional persistance length ``\\zeta``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_persistance_length::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_persistance_length,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b`` and nondimensional persistance length ``\\zeta``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_persistance_length::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_wlc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_persistance_length_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_persistance_length,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the WLC model in the isometric ensemble approximated using a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function WLC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    persistance_length::Float64,
)
    nondimensional_persistance_length = persistance_length / number_of_links / link_length
    return WLC(
        number_of_links,
        link_length,
        hinge_mass,
        persistance_length,
        nondimensional_persistance_length,
        (end_to_end_length, temperature) -> gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            persistance_length,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> gibbs_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            persistance_length,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_gibbs_free_energy(
            number_of_links,
            link_length,
            persistance_length,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_gibbs_free_energy_per_link(
            number_of_links,
            link_length,
            persistance_length,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_gibbs_free_energy(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_persistance_length,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_gibbs_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_persistance_length,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_gibbs_free_energy(
                number_of_links,
                nondimensional_persistance_length,
                nondimensional_end_to_end_length_per_link,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_gibbs_free_energy_per_link(
                number_of_links,
                nondimensional_persistance_length,
                nondimensional_end_to_end_length_per_link,
            ),
    )
end

end
