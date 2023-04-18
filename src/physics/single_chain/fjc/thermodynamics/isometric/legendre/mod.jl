"""
The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using .......Polymers: PROJECT_ROOT
using .....SingleChain: ONE, ZERO, POINTS, integrate

"""
The structure of the thermodynamics of the FJC model in the isometric ensemble approximated using a Legendre transformation.
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
The expected force as a function ``f`` of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
f(\\xi, T) \\sim \\frac{kT}{\\ell_b}\\,\\mathcal{L}^{-1}\\left(\\frac{\\xi}{N_b\\ell_b}\\right) \\quad \\text{for } N_b\\gg 1,
```

where ``\\mathcal{L}(x)=\\coth(x)-1/x`` is the Langevin function.

$(TYPEDSIGNATURES)
"""
function force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i, temperature_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
The expected nondimensional force as a function ``\\eta`` of the applied nondimensional end-to-end length per link ``\\gamma``,

```math
\\eta(\\gamma) \\sim \\mathcal{L}^{-1}(\\gamma) \\quad \\text{for } N_b\\gg 1,
```

where ``\\mathcal{L}(x)=\\coth(x)-1/x`` is the Langevin function.

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_end_to_end_length_per_link_i -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64,),
            nondimensional_end_to_end_length_per_link_i,
        ),
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The Helmholtz free energy ``\\psi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``,

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
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        end_to_end_length,
        temperature,
    )
end

"""
The Helmholtz free energy per link ``\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i, temperature_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy per link ``\\Delta\\psi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i, temperature_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional Helmholtz free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``,
given by [Buche and Silberstein](https://doi.org/10.1103/PhysRevE.102.012501) as

```math
\\vartheta(\\gamma, T) \\sim \\varphi\\left[\\mathcal{L}^{-1}(\\gamma), T\\right] + \\gamma\\mathcal{L}^{-1}(\\gamma) \\quad \\text{for } N_b\\gg 1,
```

where ``\\mathcal{L}(x)=\\coth(x)-1/x`` is the Langevin function.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64),
            number_of_links_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
given by [Buche and Silberstein](https://doi.org/10.1016/j.jmps.2021.104593) as

```math
\\Delta\\vartheta(\\gamma) \\sim \\gamma\\mathcal{L}^{-1}(\\gamma) + \\ln\\left\\{\\frac{\\mathcal{L}^{-1}(\\gamma)}{\\sinh[\\mathcal{L}^{-1}(\\gamma)]}\\right\\} \\quad \\text{for } N_b\\gg 1,
```

where ``\\mathcal{L}(x)=\\coth(x)-1/x`` is the Langevin function.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_end_to_end_length_per_link_i -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64,),
            nondimensional_end_to_end_length_per_link_i,
        ),
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The equilibrium probability density of end-to-end vectors ``P_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
P_\\mathrm{eq}(\\xi) = \\frac{e^{-\\beta\\psi(\\xi, T)}}{4\\pi\\int e^{-\\beta\\psi(\\xi', T)} \\,{\\xi'}{}^2 d\\xi'}.
```

$(TYPEDSIGNATURES)
"""
function equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_equilibrium_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
    )
end

"""
The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link ``\\mathscr{P}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``,

```math
\\mathscr{P}_\\mathrm{eq}(\\gamma) = \\frac{e^{-\\Delta\\vartheta(\\gamma)}}{4\\pi\\int e^{-\\Delta\\vartheta(\\gamma')} \\,{\\gamma'}{}^2 d\\gamma'}.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            normalization_nondimensional_equilibrium_distribution,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The equilibrium probability density of end-to-end lengths ``g_\\mathrm{eq}`` as a function of the end-to-end length ``\\xi``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
g_\\mathrm{eq}(\\xi) = 4\\pi\\xi^2 P_\\mathrm{eq}(\\xi).
```

$(TYPEDSIGNATURES)
"""
function equilibrium_radial_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_equilibrium_radial_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
    )
end

"""
The nondimensional equilibrium probability density of nondimensional end-to-end lenghts per link ``\\mathscr{g}_\\mathrm{eq}`` as a function of the nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``,

```math
\\mathscr{g}_\\mathrm{eq}(\\gamma) = 4\\pi\\gamma^2 \\mathscr{P}_\\mathrm{eq}(\\gamma).
```

$(TYPEDSIGNATURES)
"""
function nondimensional_equilibrium_radial_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    normalization_nondimensional_equilibrium_distribution::Float64,
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_equilibrium_radial_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            normalization_nondimensional_equilibrium_distribution,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The Gibbs free energy ``\\varphi`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``,

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
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        end_to_end_length,
        temperature,
    )
end

"""
The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Helmholtz free energy ``\\Delta\\varphi\\equiv\\varphi(\\xi,T)-\\varphi(0,T)`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i, temperature_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i, temperature_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64),
            number_of_links_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The nondimensional relative Helmholtz free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma`` and temperature ``T``,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_legendre_nondimensional_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64),
            number_of_links_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the isometric ensemble approximated using a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    normalization_nondimensional_equilibrium_distribution = integrate(
        nondimensional_end_to_end_length_per_link ->
            nondimensional_equilibrium_radial_distribution(
                number_of_links,
                1.0,
                nondimensional_end_to_end_length_per_link,
            ),
        ZERO,
        ONE,
        POINTS,
    )
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        normalization_nondimensional_equilibrium_distribution,
        (end_to_end_length, temperature) ->
            force(number_of_links, link_length, end_to_end_length, temperature),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_force(nondimensional_end_to_end_length_per_link),
        (end_to_end_length, temperature) -> helmholtz_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy(
            number_of_links,
            link_length,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_helmholtz_free_energy_per_link(
            number_of_links,
            link_length,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_helmholtz_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_helmholtz_free_energy(
                number_of_links,
                nondimensional_end_to_end_length_per_link,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_helmholtz_free_energy_per_link(
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length) -> equilibrium_distribution(
            number_of_links,
            link_length,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length,
        ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_distribution(
                number_of_links,
                normalization_nondimensional_equilibrium_distribution,
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length) -> equilibrium_radial_distribution(
            number_of_links,
            link_length,
            normalization_nondimensional_equilibrium_distribution,
            end_to_end_length,
        ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_radial_distribution(
                number_of_links,
                normalization_nondimensional_equilibrium_distribution,
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length, temperature) -> gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> gibbs_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_gibbs_free_energy(
            number_of_links,
            link_length,
            end_to_end_length,
            temperature,
        ),
        (end_to_end_length, temperature) -> relative_gibbs_free_energy_per_link(
            number_of_links,
            link_length,
            end_to_end_length,
            temperature,
        ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_gibbs_free_energy(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link, temperature) ->
            nondimensional_gibbs_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                nondimensional_end_to_end_length_per_link,
                temperature,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_gibbs_free_energy(
                number_of_links,
                nondimensional_end_to_end_length_per_link,
            ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_relative_gibbs_free_energy_per_link(
                number_of_links,
                nondimensional_end_to_end_length_per_link,
            ),
    )
end

end
