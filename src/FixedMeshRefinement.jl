module FixedMeshRefinement

using LinearAlgebra
using StaticArrays
using OffsetArrays

include("utils.jl")
include("grid.jl")
include("prolongation.jl")
include("restriction.jl")
include("rk.jl")
include("step.jl")

end
