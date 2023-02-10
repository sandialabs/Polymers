"""
The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach.
"""
module Asymptotic

using DocStringExtensions
using .......Polymers: PROJECT_ROOT

include("./strong_potential/mod.jl")
include("./weak_potential/mod.jl")

"""
The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach.

$(FIELDS)
"""
struct FJC
    """
    The number of links in the chain.
    """
    number_of_links::UInt8
    """
    The length of each link in the chain in units of nm.
    """
    link_length::Float64
    """
    The number of links in the chain.
    """
    hinge_mass::Float64
    """
    The thermodynamic functions of the model in the modified canonical ensemble approximated using an asymptotic approach approximated using an asymptotic approach valid for strong potentials.
    """
    strong_potential::Any
    """
    The thermodynamic functions of the model in the modified canonical ensemble approximated using an asymptotic approach approximated using an asymptotic approach valid for weak potentials.
    """
    weak_potential::Any
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        StrongPotential.FJC(number_of_links, link_length, hinge_mass),
        WeakPotential.FJC(number_of_links, link_length, hinge_mass),
    )
end

end
