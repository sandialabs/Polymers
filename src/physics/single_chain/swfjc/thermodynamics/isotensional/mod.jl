"""
The square-well freely-jointed chain (SWFJC) model thermodynamics in the isotensional ensemble.
"""
module Isotensional

using DocStringExtensions
using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

"""
The structure of the thermodynamics of the SWFJC model in the isotensional ensemble.

$(FIELDS)
"""
struct SWFJC
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
    The width of the well in units of nm.
    """
    well_width::Float64
    """
    The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    """
    legendre::Any
end

"""
Initializes and returns an instance of the thermodynamics of the SWFJC model in the isotensional ensemble.

$(TYPEDSIGNATURES)
"""
function SWFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    well_width::Float64,
)
    return SWFJC(
        number_of_links,
        link_length,
        hinge_mass,
        well_width,
        Legendre.SWFJC(number_of_links, link_length, hinge_mass, well_width),
    )
end

end
