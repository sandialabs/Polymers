"""
The Morse potential freely-jointed chain (Morse-FJC) single-chain model thermodynamics.
"""
module Thermodynamics

using DocStringExtensions

include("isometric/mod.jl")
include("isotensional/mod.jl")

"""
The structure of the thermodynamics of the Morse-FJC model.

$(FIELDS)
"""
struct MORSEFJC
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
    The stiffness of each link in the chain ``k_0`` in units of J/(mol⋅nm^2).
    """
    link_stiffness::Float64
    """
    The energy of each link in the chain ``u_0`` in units of J/mol.
    """
    link_energy::Float64
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
Initializes and returns an instance of the thermodynamics of the Morse-FJC model.

$(TYPEDSIGNATURES)
"""
function MORSEFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    link_stiffness::Float64,
    link_energy::Float64,
)
    return MORSEFJC(
        number_of_links,
        link_length,
        hinge_mass,
        link_stiffness,
        link_energy,
        Isometric.MORSEFJC(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
        ),
        Isotensional.MORSEFJC(
            number_of_links,
            link_length,
            hinge_mass,
            link_stiffness,
            link_energy,
        ),
    )
end

end
