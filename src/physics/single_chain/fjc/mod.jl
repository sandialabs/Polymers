"""
The freely-jointed chain (FJC) single-chain model.
"""
module Fjc

using DocStringExtensions
using ....Polymers: PATHSEP

include(string("thermodynamics", PATHSEP, "mod.jl"))

"""
The structure of the FJC model.

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
    The thermodynamic functions of the model.
    """
    thermodynamics::Any
end

"""
Initializes and returns an instance of the FJC model.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        Thermodynamics.FJC(number_of_links, link_length, hinge_mass),
    )
end

end
