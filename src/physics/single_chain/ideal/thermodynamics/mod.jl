"""
The ideal chain model thermodynamics.
"""
module Thermodynamics

using DocStringExtensions

include("isometric/mod.jl")
include("isotensional/mod.jl")

"""
The structure of the thermodynamics of the ideal chain model.

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
    The thermodynamic functions of the model in the isometric ensemble.
    """
    isometric::Any
    """
    The thermodynamic functions of the model in the isotensional ensemble.
    """
    isotensional::Any
end

"""
Initializes and returns an instance of the thermodynamics of the ideal chain model.

$(TYPEDSIGNATURES)
"""
function IDEAL(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return IDEAL(
        number_of_links,
        link_length,
        hinge_mass,
        Isometric.IDEAL(number_of_links, link_length, hinge_mass),
        Isotensional.IDEAL(number_of_links, link_length, hinge_mass),
    )
end

end
