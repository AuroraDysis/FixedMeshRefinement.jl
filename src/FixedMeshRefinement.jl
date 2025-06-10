module FixedMeshRefinement

using LinearAlgebra
using StaticArrays
using OffsetArrays
using FastBroadcast

include("utils.jl")
include("grid.jl")
include("prolongation.jl")
include("restriction.jl")
include("rk.jl")
include("step.jl")
include("excision.jl")

end
