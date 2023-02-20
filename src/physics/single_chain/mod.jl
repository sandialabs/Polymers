"""
Single-chain models for polymer physics.
"""
module SingleChain

using ...Polymers: PATHSEP

include("test.jl")
include(string("ideal", PATHSEP, "mod.jl"))
include(string("fjc", PATHSEP, "mod.jl"))
include(string("efjc", PATHSEP, "mod.jl"))
include(string("swfjc", PATHSEP, "mod.jl"))

ONE::Float64 = 1.0
ZERO::Float64 = 1e-6
POINTS::UInt128 = 100

function integrate(
    fun::Function,
    lower_lim::Float64,
    upper_lim::Float64,
    num_points::UInt128,
)
    dx = (upper_lim - lower_lim) / num_points
    return sum(map(fun, lower_lim .+ (0.5 .+ collect(range(0, num_points - 1))) * dx)) * dx
end

end
