"""
The freely-jointed chain (FJC) model thermodynamics.
"""
module Thermodynamics

using DocStringExtensions

include("isometric/mod.jl")
include("isotensional/mod.jl")
include("modified_canonical/mod.jl")

"""
The structure of the thermodynamics of the FJC model.

$(FIELDS)
"""
struct FJC
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
    """
    The thermodynamic functions of the model in the modified canonical ensemble.
    """
    modified_canonical::Any
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        Isometric.FJC(number_of_links, link_length, hinge_mass),
        Isotensional.FJC(number_of_links, link_length, hinge_mass),
        ModifiedCanonical.FJC(number_of_links, link_length, hinge_mass),
    )
end

end
