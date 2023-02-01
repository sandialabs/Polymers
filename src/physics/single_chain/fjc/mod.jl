module Fjc

include("./thermodynamics/mod.jl")

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    thermodynamics::Any
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.thermodynamics = Thermodynamics.FJC(number_of_links, link_length, hinge_mass)
        return fjc
    end
end

end
