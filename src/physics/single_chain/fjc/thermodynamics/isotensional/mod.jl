"""
The freely-jointed chain (FJC) model thermodynamics in the isotensional ensemble.
"""
module Isotensional

using DocStringExtensions
using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

"""
The structure of the thermodynamics of the FJC model in the isotensional ensemble.

$(FIELDS)
"""
struct FJC
    """
    The number of links in the chain.
    """
    number_of_links::UInt8
    """
    The length of each link in the chain in units of nm.
    """
    link_length::Float64
    """
    The number of links in the chain.
    """
    hinge_mass::Float64
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
end

"""
The expected end-to-end length ``\\xi`` as a function of the applied force ``f`` and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``,
given by [Rubinstein and Colby](https://global.oup.com/academic/product/polymer-physics-9780198520597) as

```math
\\xi(f, T) = -\\frac{\\partial\\varphi}{\\partial f} = N_b \\ell_b \\mathcal{L}(\\beta f\\ell_b),
```

where ``\\mathcal{L}(x)=\\coth(x)-1/x`` is the Langevin function.

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
                :physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
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
The expected end-to-end length ``\\xi`` as a function of the applied force ``f`` and temperature ``T``,
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
                :physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
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
Initializes and returns an instance of the thermodynamics of the FJC model in the isotensional ensemble.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        Legendre.FJC(number_of_links, link_length, hinge_mass),
        (force, temperature) ->
            end_to_end_length(number_of_links, link_length, force, temperature),
        (force, temperature) ->
            end_to_end_length(link_length, force, temperature),
    )
end

end
