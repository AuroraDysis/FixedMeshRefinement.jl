using BergerOligerAMR
using Test
using Aqua
using JET

@testset "BergerOligerAMR.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BergerOligerAMR)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(BergerOligerAMR; target_defined_modules = true)
    end
    # Write your tests here.
end
