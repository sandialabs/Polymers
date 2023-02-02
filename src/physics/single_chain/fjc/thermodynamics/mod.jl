module Thermodynamics

include("./isometric/mod.jl")

mutable struct FJC
    number_of_links::UInt8
    link_length::Float64
    hinge_mass::Float64
    isometric::Any
    function FJC(number_of_links::UInt8, link_length::Float64, hinge_mass::Float64)
        fjc = new(number_of_links, link_length, hinge_mass)
        fjc.isometric = Isometric.FJC(number_of_links, link_length, hinge_mass)
        return fjc
    end
end

end
