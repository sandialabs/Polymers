"""
The composite uFJC (CuFJC) single-chain model thermodynamics in the isometric ensemble approximated using a reduced asymptotic approach and a Legendre transformation.
"""
module Legendre

using DocStringExtensions

"""
The structure of the CuFJC model thermodynamics in the isometric ensemble approximated using a reduced asymptotic approach and a Legendre transformation.

$(FIELDS)
"""
struct CUFJC
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
    The number of bonds in each link.
    """
    number_of_bonds::UInt8
    """
    The stiffness of each bond in units of J/(mol⋅nm^2).
    """
    bond_stiffness::Float64
    """
    The energy of each bond in units of J/mol.
    """
    bond_energy::Float64
    """
    The scission energy of each bond in units of J/mol.
    """
    bond_scission_energy::Float64
    """
    The attempt frequency of each bond in units of 1/ns.
    """
    bond_attempt_frequency::Float64
end

"""
Initializes and returns an instance of the CuFJC model thermodynamics in the isometric ensemble approximated using a reduced asymptotic approach and a Legendre transformation.

$(TYPEDSIGNATURES)
"""
function CUFJC(
    number_of_links::UInt8,
    link_length::Float64,
    hinge_mass::Float64,
    number_of_bonds::UInt8,
    bond_stiffness::Float64,
    bond_energy::Float64,
    bond_scission_energy::Float64,
    bond_attempt_frequency::Float64,
)
    return CUFJC(
        number_of_links,
        link_length,
        hinge_mass,
        number_of_bonds,
        bond_stiffness,
        bond_energy,
        bond_scission_energy,
        bond_attempt_frequency,
    )
end

end
