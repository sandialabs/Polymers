"""
The extensible freely-jointed chain (EFJC) model thermodynamics in the isotensional ensemble approximated using an reduced asymptotic approach.
"""
module Reduced

using DocStringExtensions
using ........Polymers: PROJECT_ROOT, PATHSEP

include("./legendre/mod.jl")

"""
The structure of the thermodynamics of the EFJC model in the isotensional ensemble approximated using an reduced asymptotic approach.

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
    The thermodynamic functions of the model in the isotensional ensemble approximated using an reduced asymptotic approach.
    """
    legendre::Any
end

"""
Initializes and returns an instance of the thermodynamics of the EFJC model in the isotensional ensemble approximated using an reduced asymptotic approach.

$(TYPEDSIGNATURES)
"""
function EFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
)
    return EFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        Legendre.EFJC(number_of_links, link_length, hinge_mass, link_stiffness),
    )
end

end
