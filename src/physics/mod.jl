module Physics

include("./single_chain/mod.jl")

"""
The Boltzmann constant in units of J/(mol⋅K).
"""
BOLTZMANN_CONSTANT::Float64 = 8.314462618

"""
The Planck constant in units of J⋅ns/mol.
"""
PLANCK_CONSTANT::Float64 = 0.06350779923502961

end
