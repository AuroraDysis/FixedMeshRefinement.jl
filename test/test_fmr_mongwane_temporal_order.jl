using Test

const FMR_STANDALONE_SOURCE = normpath(joinpath(@__DIR__, "..", "src"))
const FMR_STANDALONE = isfile(joinpath(FMR_STANDALONE_SOURCE, "FixedMeshRefinement.jl"))
const FMR_SOURCE_DIR = let
    vendored = normpath(joinpath(@__DIR__, "..", "..", "src_fdm", "FixedMeshRefinement"))
    if FMR_STANDALONE
        FMR_STANDALONE_SOURCE
    elseif isfile(joinpath(vendored, "FixedMeshRefinement.jl"))
        vendored
    else
        error("FixedMeshRefinement source directory not found")
    end
end

if FMR_STANDALONE
    @eval using FixedMeshRefinement
    @eval using .FixedMeshRefinement.OffsetArrays
    @eval using .FixedMeshRefinement.FastBroadcast
else
    @eval using OffsetArrays
    @eval using FastBroadcast
end

erf(x::Float64) = ccall(:erf, Cdouble, (Cdouble,), x)
typetol(::Type{T}) where {T<:AbstractFloat} = eps(T)^(2//3)
function isapprox_tol(a::T, b::T) where {T<:AbstractFloat}
    tol = typetol(T)
    return isapprox(a, b; atol=tol, rtol=tol)
end

include(joinpath(FMR_SOURCE_DIR, "grid.jl"))
include(joinpath(FMR_SOURCE_DIR, "prolongation.jl"))
include(joinpath(FMR_SOURCE_DIR, "restriction.jl"))
include(joinpath(FMR_SOURCE_DIR, "rk.jl"))
include(joinpath(FMR_SOURCE_DIR, "step.jl"))

function run_mongwane_temporal_mms(nsteps::Int)
    t_final = 0.25
    base_points = 17
    base_dx = 1 / (base_points - 1)
    coarse_dt = t_final / nsteps
    grid = Grid(
        1,
        base_points,
        [(0.0, 1.0), (0.25, 0.75)],
        0,
        2;
        num_transition_points=0,
        spatial_interpolation_order=1,
        cfl=coarse_dt / base_dx,
    )

    for level in grid.levels
        fill!(get_state(level), 1.0)
    end

    function rhs!(level, du, u, _params, _time)
        indices = get_rk4_evaluation_indices(level)
        if level.index == 1
            @views du[:, indices] .= u[:, indices]
        else
            right_ghost = first(get_boundary_indices(level)[2])
            @views du[:, indices] .= u[1, right_ghost]
        end
        return nothing
    end

    for _ in 1:nsteps
        step!(grid, rhs!, nothing; mongwane=true)
    end

    fine = grid.levels[2]
    values = @view get_state(fine)[1, 1:(fine.num_interior_points)]
    error = maximum(abs, values .- exp(t_final))
    noise_floor = maximum(abs, values .- first(values))
    return (; coarse_dt, error, noise_floor)
end

@testset "Mongwane temporal MMS is fourth order" begin
    results = run_mongwane_temporal_mms.([4, 8, 16, 32])
    errors = getproperty.(results, :error)
    orders = log2.(errors[1:(end - 1)] ./ errors[2:end])
    @info "Mongwane temporal MMS" results orders

    @test all(result -> result.noise_floor == 0.0, results)
    @test all(order -> 3.9 <= order <= 4.05, orders)
end
