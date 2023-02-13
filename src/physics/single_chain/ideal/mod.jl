"""
The ideal single-chain model.
"""
module Ideal

using DocStringExtensions

include("./thermodynamics/mod.jl")

"""
The structure of the ideal chain model.

$(FIELDS)
"""
struct IDEAL
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
    The thermodynamic functions of the model.
    """
    thermodynamics::Any
end

"""
Initializes and returns an instance of the ideal chain model.

$(TYPEDSIGNATURES)
"""
function IDEAL(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return IDEAL(
        number_of_links,
        link_length,
        hinge_mass,
        Thermodynamics.IDEAL(number_of_links, link_length, hinge_mass),
    )
end

end
