"""
The extensible freely-jointed chain (EFJC) single-chain model.
"""
module Efjc

using DocStringExtensions

include("./thermodynamics/mod.jl")

"""
The structure of the EFJC model.

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
    The number of links in the chain.
    """
    hinge_mass::Float64
    """
    The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    """
    link_stiffness::Float64
    """
    The thermodynamic functions of the model.
    """
    thermodynamics::Any
end

"""
Initializes and returns an instance of the EFJC model.

$(TYPEDSIGNATURES)
"""
function EFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
)
    return EFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        Thermodynamics.EFJC(number_of_links, link_length, hinge_mass, link_stiffness),
    )
end

end
