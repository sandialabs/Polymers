"""
The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble.
"""
module Isometric

using DocStringExtensions
using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

"""
$(TYPEDSIGNATURES)
"""
function force(
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
                :physics_single_chain_fjc_thermodynamics_isometric_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            hinge_mass_i,
            link_length_i,
            end_to_end_length_i,
            temperature_i,
        ),
        number_of_links,
        hinge_mass,
        link_length,
        end_to_end_length,
        temperature,
    )
end

"""
$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_end_to_end_length_per_link::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_end_to_end_length_per_link_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_isometric_nondimensional_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64),
            number_of_links_i,
            hinge_mass_i,
            link_length_i,
            nondimensional_end_to_end_length_per_link_i,
        ),
        number_of_links,
        hinge_mass,
        link_length,
        nondimensional_end_to_end_length_per_link,
    )
end

"""
The structure of the thermodynamics of the FJC model in the isometric ensemble.

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
    The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    """
    legendre::Any
    """
    The expected force ``f`` as a function of the applied end-to-end length ``\\xi`` and temperature ``T``.
    """
    force::Function
    """
    The expected nondimensional force ``\\eta`` as a function of the applied nondimensional end-to-end length per link ``\\gamma``.
    """
    nondimensional_force::Function
    # helmholtz_free_energy::Function
    # helmholtz_free_energy_per_link::Function
    # relative_helmholtz_free_energy::Function
    # relative_helmholtz_free_energy_per_link::Function
    # nondimensional_helmholtz_free_energy::Function
    # nondimensional_helmholtz_free_energy_per_link::Function
    # nondimensional_relative_helmholtz_free_energy::Function
    # nondimensional_relative_helmholtz_free_energy_per_link::Function
    # equilibrium_distribution::Function
    # nondimensional_equilibrium_distribution::Function
    # equilibrium_radial_distribution::Function
    # nondimensional_equilibrium_radial_distribution::Function
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the isometric ensemble.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        Legendre.FJC(number_of_links, link_length, hinge_mass),
        (end_to_end_length, temperature) ->
            force(number_of_links, link_length, hinge_mass, end_to_end_length, temperature),
        (nondimensional_end_to_end_length_per_link) -> nondimensional_force(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_end_to_end_length_per_link,
        ),
    )
end

end
