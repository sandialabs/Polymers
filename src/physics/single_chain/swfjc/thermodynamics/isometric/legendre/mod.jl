"""
The square-well freely-jointed chain (SWFJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using Polymers_jll
using .......Polymers: PROJECT_ROOT
using .....SingleChain: ONE, ZERO, POINTS, integrate

"""
The structure of the thermodynamics of the SWFJC model in the isometric ensemble approximated using a Legendre transformation.
$(FIELDS)
"""
struct SWFJC
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
    The width of the well ``w`` in units of nm.
    """
    well_width::Float64
    normalization_nondimensional_equilibrium_distribution::Float64
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
    """
    The equilibrium probability density of end-to-end vectors ``P_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``.
    """
    equilibrium_distribution::Function
    """
    The nondimensional equilibrium probability density of end-to-end vectors ``\\mathscr{P}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_equilibrium_distribution::Function
    """
    The equilibrium probability density of end-to-end lengths ``g_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``.
    """
    equilibrium_radial_distribution::Function
    """
    The nondimensional equilibrium probability density of end-to-end lengths ``\\mathscr{g}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_equilibrium_radial_distribution::Function
end

"""
The expected force as a function ``f`` of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``.

$(TYPEDSIGNATURES)
"""
function force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_force,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        well_width,
        end_to_end_length,
        temperature,
    )
end

"""
The expected nondimensional force as a function ``\\eta`` of the applied nondimensional end-to-end length per link ``\\gamma``,
parameterized by the link length ``\\ell_b`` and well width ``w``.

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, well_width_i, nondimensional_end_to_end_length_per_link_i) ->
            ccall(
                (
                    :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_force,
                    Polymers_jll.libpolymers,
                ),
                Float64,
                (Float64, Float64, Float64),
                link_length_i,
                well_width_i,
                nondimensional_end_to_end_length_per_link_i,
            ),
        link_length,
        well_width,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The Helmholtz free energy ``\\psi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, well width ``w``, and hinge mass ``m``,

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
    well_width::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        end_to_end_length,
        temperature,
    )
end

"""
The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, well width ``w``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        well_width,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        well_width,
        end_to_end_length,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, well width ``w``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, well width ``w``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy_per_link,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        link_length,
        well_width,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the link length ``\\ell_b`` and well width ``w``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (link_length_i, well_width_i, nondimensional_end_to_end_length_per_link_i) ->
            ccall(
                (
                    :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
                    Polymers_jll.libpolymers,
                ),
                Float64,
                (Float64, Float64, Float64),
                link_length_i,
                well_width_i,
                nondimensional_end_to_end_length_per_link_i,
            ),
        link_length,
        well_width,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The equilibrium probability density of end-to-end vectors ``P_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``,

```math
P_\\mathrm{eq}(\\xi) = \\frac{e^{-\\beta\\psi(\\xi, T)}}{4\\pi\\int e^{-\\beta\\psi(\\xi', T)} \\,{\\xi'}{}^2 d\\xi'}.
```

$(TYPEDSIGNATURES)
"""
function equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, well_width_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_equilibrium_distribution,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length_i,
        ),
        number_of_links,
        link_length,
        well_width,
        end_to_end_length,
    )
end

"""
The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link ``\\mathscr{P}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``,

```math
\\mathscr{P}_\\mathrm{eq}(\\gamma) = \\frac{e^{-\\Delta\\vartheta(\\gamma)}}{4\\pi\\int e^{-\\Delta\\vartheta(\\gamma')} \\,{\\gamma'}{}^2 d\\gamma'}.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_distribution,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            normalization_nondimensional_equilibrium_distribution,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        link_length,
        well_width,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The equilibrium probability density of end-to-end lengths ``g_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``,

```math
g_\\mathrm{eq}(\\xi) = 4\\pi\\xi^2 P_\\mathrm{eq}(\\xi).
```

$(TYPEDSIGNATURES)
"""
function equilibrium_radial_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, well_width_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_equilibrium_radial_distribution,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length_i,
        ),
        number_of_links,
        link_length,
        well_width,
        end_to_end_length,
    )
end

"""
The nondimensional equilibrium probability density of nondimensional end-to-end lenghts per link ``\\mathscr{g}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and well width ``w``,

```math
\\mathscr{g}_\\mathrm{eq}(\\gamma) = 4\\pi\\gamma^2 \\mathscr{P}_\\mathrm{eq}(\\gamma).
```

$(TYPEDSIGNATURES)
"""
function nondimensional_equilibrium_radial_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    well_width::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            well_width_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_swfjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_radial_distribution,
                Polymers_jll.libpolymers,
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            well_width_i,
            normalization_nondimensional_equilibrium_distribution,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        link_length,
        well_width,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the SWFJC model in the isometric ensemble approximated using a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function SWFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    well_width::Float64,
)
    normalization_nondimensional_equilibrium_distribution = integrate(
        nondimensional_end_to_end_length_per_link ->
            nondimensional_equilibrium_radial_distribution(
                number_of_links,
                link_length,
                well_width,
                1.0,
                nondimensional_end_to_end_length_per_link,
            ),
        ZERO,
        ONE * (1.0 + well_width / link_length),
        POINTS,
    )
    return SWFJC(
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        normalization_nondimensional_equilibrium_distribution,
        (end_to_end_length, temperature) ->
            force(number_of_links, link_length, well_width, end_to_end_length, temperature),
        (nondimensional_end_to_end_length_per_link) -> nondimensional_force(
            link_length,
            well_width,
            nondimensional_end_to_end_length_per_link,
        ),
        (end_to_end_length, temperature) -> helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            well_width,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            well_width,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy(
            number_of_links,
            link_length,
            well_width,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            well_width,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy(
                number_of_links,
                link_length,
                hinge_mass,
                well_width,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                well_width,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_helmholtz_free_energy(
                number_of_links,
                link_length,
                well_width,
                nondimensional_end_to_end_length_per_link,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_helmholtz_free_energy_per_link(
                link_length,
                well_width,
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length) -> equilibrium_distribution(
            number_of_links,
            link_length,
            well_width,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length,
        ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_distribution(
                number_of_links,
                link_length,
                well_width,
                normalization_nondimensional_equilibrium_distribution,
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length) -> equilibrium_radial_distribution(
            number_of_links,
            link_length,
            well_width,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length,
        ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_radial_distribution(
                number_of_links,
                link_length,
                well_width,
                normalization_nondimensional_equilibrium_distribution,
                nondimensional_end_to_end_length_per_link,
            ),
    )
end

end
