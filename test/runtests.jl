using Polymers
using Test

@testset "Polymers.jl" begin
    @test Polymers.get_name() != "Hello world!"
end
