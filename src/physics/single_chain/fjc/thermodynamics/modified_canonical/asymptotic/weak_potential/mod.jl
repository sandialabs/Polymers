"""
The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.
"""
module WeakPotential

using DocStringExtensions
using ........Polymers: PROJECT_ROOT

"""
The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.

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
    The expected end-to-end length ``\\xi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    end_to_end_length::Function
    """
    The expected end-to-end length per link ``\\xi/N_b=\\ell_b\\gamma`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    end_to_end_length_per_link::Function
    """
    The expected nondimensional end-to-end length ``N_b\\gamma=\\xi/\\ell_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_end_to_end_length::Function
    """
    The expected nondimensional end-to-end length per link ``\\gamma\\equiv\\xi/N_b\\ell_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_end_to_end_length_per_link::Function
    """
    The expected force ``f`` as a function of the applied potential distance and potential stiffness.
    """
    force::Function
    """
    The expected nondimensional force ``\\eta`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_force::Function
    """
    The Gibbs free energy ``\\varphi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    gibbs_free_energy::Function
    """
    The Gibbs free energy per link ``\\varphi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    gibbs_free_energy_per_link::Function
    """
    The relative Gibbs free energy ``\\Delta\\varphi\\equiv\\varphi(f,T)-\\varphi(0,T)`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    relative_gibbs_free_energy::Function
    """
    The relative Gibbs free energy per link ``\\Delta\\varphi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``.
    """
    relative_gibbs_free_energy_per_link::Function
    """
    The nondimensional Gibbs free energy ``N_b\\varrho=\\beta\\varphi`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``.
    """
    nondimensional_gibbs_free_energy::Function
    """
    The nondimensional Gibbs free energy per link ``\\varrho\\equiv\\beta\\varphi/N_b`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``.
    """
    nondimensional_gibbs_free_energy_per_link::Function
    """
    The nondimensional relative Gibbs free energy ``N_b\\Delta\\varrho=\\beta\\Delta\\varphi`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_relative_gibbs_free_energy::Function
    """
    The nondimensional relative Gibbs free energy per link ``\\Delta\\varrho\\equiv\\beta\\Delta\\varphi/N_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    """
    nondimensional_relative_gibbs_free_energy_per_link::Function
end

"""
The expected end-to-end length ``\\xi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_end_to_end_length,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The expected end-to-end length per link ``\\xi/N_b=\\ell_b\\gamma`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function end_to_end_length_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_end_to_end_length_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The expected nondimensional end-to-end length ``N_b\\gamma\\equiv\\xi/\\ell_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_end_to_end_length,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The expected nondimensional end-to-end length ``\\gamma\\equiv\\xi/N_b\\ell_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_end_to_end_length_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_end_to_end_length_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The expected force ``f`` as a function of the applied potential distance and potential stiffness.

$(TYPEDSIGNATURES)
"""
function force(
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (potential_distance_i, potential_stiffness_i) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (Float64, Float64),
            potential_distance_i,
            potential_stiffness_i,
        ),
        potential_distance,
        potential_stiffness,
    )
end

"""
The expected nondimensional force ``\\eta`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness.

$(TYPEDSIGNATURES)
"""
function nondimensional_force(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_force,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The gibbs free energy ``\\psi`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The gibbs free energy per link ``\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The relative gibbs free energy ``\\Delta\\psi\\equiv\\psi(\\xi,T)-\\psi(0,T)`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_relative_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The relative gibbs free energy per link ``\\Delta\\psi/N_b`` as a function of the applied potential distance, potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b`` and link length ``\\ell_b``.

$(TYPEDSIGNATURES)
"""
function relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    potential_distance::Union{Float64,Vector,Matrix,Array},
    potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            potential_distance_i,
            potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        potential_distance,
        potential_stiffness,
        temperature,
    )
end

"""
The nondimensional gibbs free energy ``N_b\\vartheta=\\beta\\psi`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
        temperature,
    )
end

"""
The nondimensional gibbs free energy per link ``\\vartheta\\equiv\\beta\\psi/N_b`` as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature ``T``,
parameterized by the number of links ``N_b``, link length ``\\ell_b``, and hinge mass ``m``.

$(TYPEDSIGNATURES)
"""
function nondimensional_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    link_length::Union{Float64,Vector,Matrix,Array},
    hinge_mass::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
    temperature::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64, Float64, Float64, Float64),
            number_of_links_i,
            link_length_i,
            hinge_mass_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
            temperature_i,
        ),
        number_of_links,
        link_length,
        hinge_mass,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
        temperature,
    )
end

"""
The nondimensional relative gibbs free energy ``N_b\\Delta\\vartheta=\\beta\\Delta\\psi`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_relative_gibbs_free_energy,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
The nondimensional relative gibbs free energy per link ``\\Delta\\vartheta\\equiv\\beta\\Delta\\psi/N_b`` as a function of the applied nondimensional potential distance and nondimensional potential stiffness,
parameterized by the number of links ``N_b``.

$(TYPEDSIGNATURES)
"""
function nondimensional_relative_gibbs_free_energy_per_link(
    number_of_links::Union{UInt8,Vector,Matrix,Array},
    nondimensional_potential_distance::Union{Float64,Vector,Matrix,Array},
    nondimensional_potential_stiffness::Union{Float64,Vector,Matrix,Array},
)::Union{Float64,Vector,Matrix,Array}
    return broadcast(
        (
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ) -> ccall(
            (
                :physics_single_chain_fjc_thermodynamics_modified_canonical_asymptotic_weak_potential_nondimensional_relative_gibbs_free_energy_per_link,
                string(PROJECT_ROOT, "target/debug/libpolymers"),
            ),
            Float64,
            (UInt8, Float64, Float64),
            number_of_links_i,
            nondimensional_potential_distance_i,
            nondimensional_potential_stiffness_i,
        ),
        number_of_links,
        nondimensional_potential_distance,
        nondimensional_potential_stiffness,
    )
end

"""
Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.

$(TYPEDSIGNATURES)
"""
function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
    return FJC(
        number_of_links,
        link_length,
        hinge_mass,
        (potential_distance, potential_stiffness, temperature) -> end_to_end_length(
            number_of_links,
            link_length,
            potential_distance,
            potential_stiffness,
            temperature,
        ),
        (potential_distance, potential_stiffness, temperature) ->
            end_to_end_length_per_link(
                number_of_links,
                link_length,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_end_to_end_length(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_end_to_end_length_per_link(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (potential_distance, potential_stiffness) ->
            force(potential_distance, potential_stiffness),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_force(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (potential_distance, potential_stiffness, temperature) -> gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            potential_distance,
            potential_stiffness,
            temperature,
        ),
        (potential_distance, potential_stiffness, temperature) ->
            gibbs_free_energy_per_link(
                number_of_links,
                link_length,
                hinge_mass,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (potential_distance, potential_stiffness, temperature) ->
            relative_gibbs_free_energy(
                number_of_links,
                link_length,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (potential_distance, potential_stiffness, temperature) ->
            relative_gibbs_free_energy_per_link(
                number_of_links,
                link_length,
                potential_distance,
                potential_stiffness,
                temperature,
            ),
        (
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ) -> nondimensional_gibbs_free_energy(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ),
        (
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ) -> nondimensional_gibbs_free_energy_per_link(
            number_of_links,
            link_length,
            hinge_mass,
            nondimensional_potential_distance,
            nondimensional_potential_stiffness,
            temperature,
        ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_relative_gibbs_free_energy(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
        (nondimensional_potential_distance, nondimensional_potential_stiffness) ->
            nondimensional_relative_gibbs_free_energy_per_link(
                number_of_links,
                nondimensional_potential_distance,
                nondimensional_potential_stiffness,
            ),
    )
end

end
