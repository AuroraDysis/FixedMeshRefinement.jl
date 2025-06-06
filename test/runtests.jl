using FixedMeshRefinement
using Test
using Aqua
using JET

@testset "FixedMeshRefinement.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(FixedMeshRefinement)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(FixedMeshRefinement; target_defined_modules = true)
    end
    # Write your tests here.
end
