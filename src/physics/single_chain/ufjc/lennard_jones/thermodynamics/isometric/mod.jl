"""
The Lennard-Jones potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isometric ensemble.
"""
module Isometric

using DocStringExtensions
using .......Polymers: PROJECT_ROOT

include("asymptotic/mod.jl")

"""
The structure of the thermodynamics of the Lennard-Jones-FJC model in the isometric ensemble.

$(FIELDS)
"""
struct LENNARDJONESFJC
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
    The thermodynamic functions of the model in the isometric ensemble approximated using an asymptotic approach.
    """
    asymptotic::Any
end

"""
Initializes and returns an instance of the thermodynamics of the Lennard-Jones-FJC model in the isometric ensemble.

$(TYPEDSIGNATURES)
"""
function LENNARDJONESFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
)
    return LENNARDJONESFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        Asymptotic.LENNARDJONESFJC(number_of_links, link_length, hinge_mass, link_stiffness),
    )
end

end
