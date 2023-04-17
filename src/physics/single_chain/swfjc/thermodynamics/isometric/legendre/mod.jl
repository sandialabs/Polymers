"""
The square-well freely-jointed chain (SWFJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
"""
module Legendre

using DocStringExtensions
using .......Polymers: PROJECT_ROOT

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
    return SWFJC(number_of_links, link_length, hinge_mass, well_width)
end

end
