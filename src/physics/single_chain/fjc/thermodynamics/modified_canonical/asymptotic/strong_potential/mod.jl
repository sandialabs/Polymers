"""
The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
"""
module StrongPotential

using DocStringExtensions
using ........Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.

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
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(number_of_links, link_length, hinge_mass)
end

end
