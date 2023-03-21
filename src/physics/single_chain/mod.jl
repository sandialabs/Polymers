"""
Single-chain models for polymer physics.
"""
module SingleChain

include("ideal/mod.jl")
include("fjc/mod.jl")
include("efjc/mod.jl")
include("swfjc/mod.jl")

const ONE::Float64 = 1.0
const ZERO::Float64 = 1e-6
const POINTS::UInt128 = 64

function integrate(
    fun::Function,
    lower_lim::Float64,
    upper_lim::Float64,
    num_points::UInt128,
)
    dx = (upper_lim - lower_lim) / num_points
    return sum(map(fun, lower_lim .+ (0.5 .+ collect(range(0, num_points - 1))) * dx)) * dx
end

struct Parameters
    abs_tol::Float64
    rel_tol::Float64
    rel_tol_thermodynamic_limit::Float64
    log_log_tol::Float64
    log_log_scale::Float64
    number_of_loops::UInt32
    hinge_mass_reference::Float64
    hinge_mass_scale::Float64
    link_length_reference::Float64
    link_length_scale::Float64
    number_of_links_minimum::UInt8
    number_of_links_maximum::UInt8
    link_stiffness_reference::Float64
    link_stiffness_scale::Float64
    nondimensional_link_stiffness_large::Float64
    nondimensional_link_stiffness_medium::Float64
    well_width_reference::Float64
    well_width_scale::Float64
    nondimensional_well_width_small::Float64
    nondimensional_end_to_end_length_per_link_reference::Float64
    nondimensional_end_to_end_length_per_link_scale::Float64
    nondimensional_end_to_end_length_per_link_small::Float64
    nondimensional_force_reference::Float64
    nondimensional_force_scale::Float64
    nondimensional_force_small::Float64
    nondimensional_potential_distance_reference::Float64
    nondimensional_potential_distance_scale::Float64
    nondimensional_potential_distance_small::Float64
    nondimensional_potential_distance_large_1::Float64
    nondimensional_potential_distance_large_2::Float64
    nondimensional_potential_stiffness_reference::Float64
    nondimensional_potential_stiffness_scale::Float64
    nondimensional_potential_stiffness_small::Float64
    nondimensional_potential_stiffness_large::Float64
    temperature_reference::Float64
    temperature_scale::Float64
end

parameters = Parameters(
    1e-7,
    1e-5,
    1e-1,
    5e-1,
    12e-1,
    8,
    1e0,
    1e0,
    1e0,
    1e0,
    0x08,
    0x19,
    5e5,
    99e4,
    1e4,
    1e1,
    99e-2,
    5e-1,
    1e-2,
    5e-1,
    99e-2,
    25e-2,
    5e1,
    1e2,
    75e-2,
    1e0,
    2e0,
    33e-2,
    1e1,
    1e1 + 25e-1,
    5e0,
    1e1,
    1e-2,
    1e2,
    3e2,
    1e2,
)

end
