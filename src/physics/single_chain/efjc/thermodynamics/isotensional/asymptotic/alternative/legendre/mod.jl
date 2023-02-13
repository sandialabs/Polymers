"""
The extensible freely-jointed chain (EFJC) model thermodynamics in the isotensional ensemble approximated using a alternative asymptotic approach and a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using .........Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the EFJC model in the isotensional ensemble approximated using a alternative asymptotic approach and a Legendre transformation.
$(FIELDS)
"""
struct EFJC
    """
    The number of links in the chain.
    """
    number_of_links::UInt8
    """
    The length of each link in the chain in units of nm.
    """
    link_length::Float64
    """
    The mass of each hinge in the chain in units of kg/mol.
    """
    hinge_mass::Float64
    """
    The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    """
    link_stiffness::Float64
end

"""
Initializes and returns an instance of the thermodynamics of the EFJC model in the isotensional ensemble approximated using a alternative asymptotic approach and a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function EFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
)
    return EFJC(number_of_links, link_length, hinge_mass, link_stiffness)
end

end
