"""
The ideal chain model thermodynamics in the isometric ensemble.
"""
module Isometric

using DocStringExtensions
using ......Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the ideal chain model in the isometric ensemble.

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
    The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``.
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
The expected force ``f`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
f(\\xi, T) = \\frac{\\partial\\psi}{\\partial\\xi} = \\frac{3kT\\xi}{N_b\\ell_b^2}.
```

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
                :physics_single_chain_ideal_thermodynamics_isometric_force,
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
The expected nondimensional force ``\\eta`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``,

```math
\\eta(\\gamma) = \\frac{\\partial\\vartheta}{\\partial\\gamma} = 3\\gamma.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_end_to_end_length_per_link_i -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_force,
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
\\psi(\\xi, T) = -kT\\ln Q(\\xi, T).
```

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
                :physics_single_chain_ideal_thermodynamics_isometric_helmholtz_free_energy,
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
                :physics_single_chain_ideal_thermodynamics_isometric_helmholtz_free_energy_per_link,
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
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,

```math
\\Delta\\psi(\\xi, T) = kT\\ln\\left[\\frac{P_\\mathrm{eq}(0)}{P_\\mathrm{eq}(\\xi)}\\right] = \\frac{3kT\\xi^2}{2N_b\\ell_b^2}.
```

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
                :physics_single_chain_ideal_thermodynamics_isometric_relative_helmholtz_free_energy,
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
                :physics_single_chain_ideal_thermodynamics_isometric_relative_helmholtz_free_energy_per_link,
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
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_helmholtz_free_energy,
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
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_helmholtz_free_energy_per_link(
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_helmholtz_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64, Float64, Float64),
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
            temperature_i,
        ),
        link_length,
        hinge_mass,
        nondimensional_end_to_end_length_per_link,
        temperature,
    )
end

"""
The nondimensional relative Helmholtz free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``,

```math
\\beta\\Delta\\psi(\\gamma) = \\ln\\left[\\frac{\\mathscr{P}_\\mathrm{eq}(0)}{\\mathscr{P}_\\mathrm{eq}(\\gamma)}\\right] = \\frac{3}{2}\\,N_b\\gamma^2.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy,
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
The nondimensional relative Helmholtz free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``,
parameterized by the number of links ``N_b``,

```math
\\Delta\\vartheta(\\gamma) = \\ln\\left[\\frac{\\mathscr{P}_\\mathrm{eq}(0)}{\\mathscr{P}_\\mathrm{eq}(\\gamma)}\\right]^{1/N_b} = \\frac{3}{2}\\,\\gamma^2.
```

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_helmholtz_free_energy_per_link(
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        nondimensional_end_to_end_length_per_link_i -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy_per_link,
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
P_\\mathrm{eq}(\\xi) = \\frac{e^{-\\beta\\psi(\\xi, T)}}{4\\pi\\int e^{-\\beta\\psi(\\xi', T)} \\,{\\xi'}{}^2 d\\xi'} = \\left(\\frac{3}{2\\pi N_b\\ell_b^2}\\right)^{3/2}\\exp\\left(-\\frac{3\\xi^2}{2N_b\\ell_b^2}\\right).
```

$(TYPEDSIGNATURES)
"""
function equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_equilibrium_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            link_length_i,
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
\\mathscr{P}_\\mathrm{eq}(\\gamma) = \\frac{e^{-\\Delta\\vartheta(\\gamma)}}{4\\pi\\int e^{-\\Delta\\vartheta(\\gamma')} \\,{\\gamma'}{}^2 d\\gamma'} = \\left(\\frac{3}{2\\pi N_b}\\right)^{3/2}\\exp\\left(-\\frac{3}{2}\\,N_b\\gamma^2\\right).
```

$(TYPEDSIGNATURES)
"""
function nondimensional_equilibrium_distribution(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_equilibrium_distribution,
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
    end_to_end_length::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, link_length_i, end_to_end_length_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_equilibrium_radial_distribution,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            link_length_i,
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
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (number_of_links_i, nondimensional_end_to_end_length_per_link_i) -> ccall(
            (
                :physics_single_chain_ideal_thermodynamics_isometric_nondimensional_equilibrium_radial_distribution,
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
Initializes and returns an instance of the thermodynamics of the ideal chain model in the isometric ensemble.

$(TYPEDSIGNATURES)
"""
function IDEAL(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return IDEAL(
        number_of_links,
        link_length,
        hinge_mass,
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
        (end_to_end_length) ->
            equilibrium_distribution(number_of_links, link_length, end_to_end_length),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_distribution(
                number_of_links,
                nondimensional_end_to_end_length_per_link,
            ),
        (end_to_end_length) -> equilibrium_radial_distribution(
            number_of_links,
            link_length,
            end_to_end_length,
        ),
        (nondimensional_end_to_end_length_per_link) ->
            nondimensional_equilibrium_radial_distribution(
                number_of_links,
                nondimensional_end_to_end_length_per_link,
            ),
    )
end

end
