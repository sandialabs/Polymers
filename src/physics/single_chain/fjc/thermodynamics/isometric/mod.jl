module Isometric

using ......Polymers: PROJECT_ROOT

include("./legendre/mod.jl")

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    legendre
    force
    nondimensional_force
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.legendre = Legendre.FJC(number_of_links, link_length, hinge_mass)
        fjc.force = (end_to_end_length, temperature) -> ccall((:fjc_thermodynamics_isometric_force, string(PROJECT_ROOT, "target/debug/libpolymers")), Float64, (UInt8, Float64, Float64, Float64, Float64), number_of_links, hinge_mass, link_length, end_to_end_length, temperature)
        fjc.nondimensional_force = nondimensional_end_to_end_length_per_link -> ccall((:fjc_thermodynamics_isometric_nondimensional_force, string(PROJECT_ROOT, "target/debug/libpolymers")), Float64, (UInt8, Float64, Float64, Float64), number_of_links, hinge_mass, link_length, nondimensional_end_to_end_length_per_link)
        return fjc
    end
end

end