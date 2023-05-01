"""
The worm-like chain (WLC) single-chain model.
"""
module Wlc

using DocStringExtensions

include("thermodynamics/mod.jl")

"""
The structure of the WLC model.

$(FIELDS)
"""
struct WLC
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
    The persistance length of the chain in units of nm.
    """
    persistance_length::Float64
    """
    The thermodynamic functions of the model.
    """
    thermodynamics::Any
end

"""
Initializes and returns an instance of the WLC model.

$(TYPEDSIGNATURES)
"""
function WLC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    persistance_length::Float64,
)
    return WLC(
        number_of_links,
        link_length,
        hinge_mass,
        persistance_length,
        Thermodynamics.WLC(number_of_links, link_length, hinge_mass, persistance_length),
    )
end

end
