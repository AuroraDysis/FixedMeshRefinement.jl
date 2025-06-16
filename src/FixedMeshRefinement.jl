module FixedMeshRefinement

using LinearAlgebra
using OffsetArrays
using BlockArrays
using FastBroadcast

include("utils.jl")
include("grid.jl")
include("prolongation.jl")
include("restriction.jl")
include("rk.jl")
include("step.jl")
include("excision.jl")
include("merge.jl")

end
