"""
The square-well freely-jointed chain (SWFJC) model thermodynamics.
"""
module Thermodynamics

using DocStringExtensions

include("./isotensional/mod.jl")

"""
The structure of the thermodynamics of the SWFJC model.

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
    The thermodynamic functions of the model in the isotensional ensemble.
    """
    isotensional::Any
end

"""
Initializes and returns an instance of the thermodynamics of the SWFJC model.

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
        Isotensional.SWFJC(number_of_links, link_length, hinge_mass, well_width),
    )
end

end
