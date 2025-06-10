export Level, Grid, shift_physical_boundary!

# Find the index of the nearest value in a sorted array
"""
    searchsortednearest(a::AbstractVector, x)

Find the index of the element in a sorted array `a` that is nearest to `x`.
"""
function searchsortednearest(a::AbstractVector, x)
    idx = searchsortedfirst(a, x)
    if idx == 1
        return 1
    elseif idx > length(a)
        return length(a)
    else
        return abs(a[idx] - x) < abs(a[idx - 1] - x) ? idx : idx - 1
    end
end

"""
    Level{NumState, NumDiagnostic}

A struct representing a single refinement level in the mesh.

# Type Parameters
- `NumState`: Number of state variables.
- `NumDiagnostic`: Number of diagnostic variables.

# Fields
- `num_interior_points::Int`: Number of interior grid points.
- `num_ghost_points::Int`: Number of ghost points on each side.
- `num_transition_points::NTuple{2,Int}`: Number of points in the transition zone for mesh refinement.
- `num_total_points::Int`: Total number of grid points (interior + ghost + buffer).
- `num_additional_points::NTuple{2,Int}`: Number of additional points on each side (ghost or buffer).
- `time_interpolation_order::Int`: Order of time interpolation.
- `spatial_interpolation_order::Int`: Order of spatial interpolation.
- `domain_box::Tuple{Float64,Float64}`: The computational domain (interior) of this level.
- `physical_domain_box::Tuple{Float64,Float64}`: The physical domain of the entire grid.
- `is_physical_boundary::NTuple{2,Bool}`: Indicates if the level boundaries are physical boundaries.
- `dx::Float64`: Grid spacing.
- `dt::Float64`: Time step.
- `t::Float64`: Current time of this level.
- `is_base_level::Bool`: True if this is the coarsest level.
- `parent_indices::UnitRange{Int}`: Range of indices in the parent level that this level covers.
- `additional_points_indices::NTuple{2,StepRange{Int,Int}}`: Indices for additional points on each side.
- `x::OffsetVector{Float64,...}`: Grid point coordinates.
- `state::Vector{OffsetMatrix{Float64,...}}`: State variables at different time levels.
- `rhs::OffsetMatrix{Float64,...}`: Right-hand side of the evolution equations.
- `tmp::OffsetMatrix{Float64,...}`: Temporary storage array.
- `k::Vector{OffsetMatrix{Float64,...}}`: Intermediate states for Runge-Kutta time integration.
- `Yn_buffer::Vector{Array{Float64,3}}`: Buffer for subcycling.
- `diag_state::OffsetMatrix{Float64,...}`: Diagnostic variables.
"""
mutable struct Level{NumState,NumDiagnostic}
    num_interior_points::Int  # num of interior grid points
    const num_ghost_points::Int  # num of ghost points on each side
    const num_transition_points::NTuple{2,Int}  # num of transition zone points
    num_total_points::Int  # num of all grid points
    const num_additional_points::NTuple{2,Int} # num of additional points on each side
    const time_interpolation_order::Int  # interpolation order in time
    const spatial_interpolation_order::Int  # interpolation order in space
    domain_box::Tuple{Float64,Float64}  # computational domain (interior)
    physical_domain_box::Tuple{Float64,Float64}  # physical domain
    const is_physical_boundary::NTuple{2,Bool}  # whether the ghost points require physical boundary condition
    const dx::Float64
    dt::Float64
    t::Float64
    const is_base_level::Bool
    parent_indices::UnitRange{Int}
    additional_points_indices::NTuple{2,StepRange{Int,Int}}

    # data
    x::OffsetVector{Float64,LinRange{Float64,Int64}}
    state::Vector{OffsetMatrix{Float64,Matrix{Float64}}} # state vectors at different time levels
    rhs::OffsetMatrix{Float64,Matrix{Float64}} # rhs of state vectors
    tmp::OffsetMatrix{Float64,Matrix{Float64}}

    # intermediate state vectors for new subcycling
    k::Vector{OffsetMatrix{Float64,Matrix{Float64}}}
    Yn_buffer::Vector{Array{Float64,3}}

    # diagnostic variables
    diag_state::OffsetMatrix{Float64,Matrix{Float64}} # state vectors for diagnostic variables

    function Level{NumState,NumDiagnostic}(
        num_interior_points,
        num_ghost_points,
        num_buffer_points,
        num_transition_points,
        time_interpolation_order,
        spatial_interpolation_order,
        domain_box,
        physical_domain_box,
        dt,
        t,
        is_base_level,
        parent_indices,
    ) where {NumState,NumDiagnostic}
        is_physical_boundary = (
            isapprox(domain_box[1], physical_domain_box[1]; rtol=typetol(Float64)),
            isapprox(domain_box[2], physical_domain_box[2]; rtol=typetol(Float64)),
        )
        num_left_additional_points =
            is_physical_boundary[1] ? num_ghost_points : num_buffer_points
        num_right_additional_points =
            is_physical_boundary[2] ? num_ghost_points : num_buffer_points
        num_total_points =
            num_interior_points + num_left_additional_points + num_right_additional_points
        dx = (domain_box[2] - domain_box[1]) / (num_interior_points - 1)
        x_min = domain_box[1] - num_left_additional_points * dx
        x_max = domain_box[2] + num_right_additional_points * dx
        x = LinRange(x_min, x_max, num_total_points)
        time_levels = max(time_interpolation_order + 1, 2)
        state = [fill(NaN, num_total_points, NumState) for _ in 1:time_levels]
        rhs = fill(NaN, num_total_points, NumState)
        tmp = fill(NaN, num_total_points, NumState)
        k = Vector{Matrix{Float64}}(undef, 4)
        Yn_buffer = Vector{Array{Float64,3}}(undef, 4)
        diag_state = fill(NaN, num_total_points, NumDiagnostic)
        for j in 1:4
            k[j] = fill(NaN, num_total_points, NumState)
            Yn_buffer[j] = fill(NaN, num_buffer_points, NumState, 2)
        end
        additional_points_indices = (
            0:-1:(-num_left_additional_points + 1),
            (num_interior_points + 1):(num_interior_points + num_right_additional_points),
        )
        offset_indices =
            (-num_left_additional_points + 1):(num_interior_points + num_right_additional_points)
        num_left_transition_points = is_physical_boundary[1] ? num_transition_points : 0
        num_right_transition_points = is_physical_boundary[2] ? num_transition_points : 0
        return new{NumState,NumDiagnostic}(
            num_interior_points,
            num_ghost_points,
            (num_left_transition_points, num_right_transition_points),
            num_total_points,
            (num_left_additional_points, num_right_additional_points),
            time_interpolation_order,
            spatial_interpolation_order,
            domain_box,
            physical_domain_box,
            is_physical_boundary,
            dx,
            dt,
            t,
            is_base_level,
            parent_indices,
            additional_points_indices,
            # data
            OffsetArray(x, offset_indices),
            [OffsetArray(state[i], offset_indices, :) for i in 1:time_levels],
            OffsetArray(rhs, offset_indices, :),
            OffsetArray(tmp, offset_indices, :),
            [OffsetArray(k[i], offset_indices, :) for i in 1:4],
            Yn_buffer,
            OffsetArray(diag_state, offset_indices, :),
        )
    end
end

"""
    cycle_state!(level::Level)

Shift the time levels of the state variables in a `Level`. `state[i]` becomes
`state[i+1]`. This is used to advance the solution in time.
"""
function cycle_state!(level::Level)
    (; state, time_interpolation_order) = level
    # shift state vectors
    for i in 1:time_interpolation_order
        state[i] .= state[i + 1]
    end
    return nothing
end

"""
    fidx2cidx(fine_level::Level, fidx::Int) -> Int

Convert a fine grid index `fidx` to the corresponding coarse grid index.
This function is not defined for the base level. An error is thrown if the fine
grid point does not align with any coarse grid point.
"""
function fidx2cidx(fine_level::Level, fidx::Int)
    (; is_base_level, parent_indices) = fine_level

    if is_base_level
        error("fidx2cidx is not defined for base level")
    end

    parent_idx_left = parent_indices[1]
    offset = fidx - 1

    if mod(offset, 2) != 0
        println("parent_indices = ", parent_indices)
        println("fidx = $fidx, offset = $offset")
        error("fidx = $fidx is not aligned with any point in parent level")
    end

    half_offset = div(offset, 2)
    return parent_idx_left + half_offset
end

"""
    Grid{NumState, NumDiagnostic}

A struct representing the entire FMR grid, which consists of multiple `Level`s.

# Type Parameters
- `NumState`: Number of state variables.
- `NumDiagnostic`: Number of diagnostic variables.

# Fields
- `num_levels::Int`: The total number of refinement levels.
- `levels::Vector{Level{NumState,NumDiagnostic}}`: A vector of `Level` objects.
- `base_dt::Float64`: The time step of the coarsest level.
- `t::Float64`: The current time of the simulation.
- `subcycling::Bool`: A boolean indicating whether subcycling in time is enabled.
"""
mutable struct Grid{NumState,NumDiagnostic}
    num_levels::Int
    levels::Vector{Level{NumState,NumDiagnostic}}
    base_dt::Float64
    t::Float64
    subcycling::Bool  # turn on subcycling or not

    function Grid{NumState,NumDiagnostic}(
        base_level_num_points,  # num of interior grid points at base level
        domain_boxes::Vector{Tuple{Float64,Float64}},
        num_ghost_points,
        num_buffer_points;
        num_transition_points=3,
        time_interpolation_order=2,
        spatial_interpolation_order=5,
        cfl=0.25,
        initial_time=0.0,
        subcycling=true,
    ) where {NumState,NumDiagnostic}
        num_levels = length(domain_boxes)
        physical_domain_box = domain_boxes[1]

        # build the coarsest level (base level)
        base_dx = (domain_boxes[1][2] - domain_boxes[1][1]) / (base_level_num_points - 1)
        base_dt = if subcycling
            cfl * base_dx
        else
            cfl * base_dx / 2^(num_levels - 1)
        end
        base_level = Level{NumState,NumDiagnostic}(
            base_level_num_points,
            num_ghost_points,
            num_buffer_points,
            num_transition_points,
            time_interpolation_order,
            spatial_interpolation_order,
            domain_boxes[1],
            physical_domain_box,
            base_dt,
            initial_time,
            true,
            (0:0),
        )
        levels = [base_level]

        # build the finer levels
        for l in 2:num_levels
            level_dx = base_dx / 2^(l - 1)
            level_dt = (subcycling ? cfl * level_dx : base_dt)
            parent_level = levels[l - 1]  # level lower than the current level (parent level)
            # if we refine parent level everywhere
            parent_level_grid_points = LinRange(
                parent_level.domain_box[1],
                parent_level.domain_box[2],
                (parent_level.num_interior_points - 1) * 2 + 1,
            )
            level_domain_box = domain_boxes[l]
            # find those two which are closest to current level boundaries using binary search
            level_min_idx = searchsortednearest(
                parent_level_grid_points, level_domain_box[1]
            )
            level_max_idx = searchsortednearest(
                parent_level_grid_points, level_domain_box[2]
            )
            level_domain = (
                parent_level_grid_points[level_min_idx],
                parent_level_grid_points[level_max_idx],
            )
            level_num_interior_points = (level_max_idx - level_min_idx) + 1
            # maps between two levels
            parent_idx_left = div(level_min_idx + 1, 2)
            parent_idx_right = div(level_max_idx + 1, 2)
            parent_indices = parent_idx_left:parent_idx_right

            # check x are aligned
            if !(
                isapprox(parent_level.x[parent_indices[1]], level_domain[1]; rtol=1e-12) &&
                isapprox(parent_level.x[parent_indices[end]], level_domain[2]; rtol=1e-12)
            )
                println("parent_indices = ", parent_indices)
                println("parent_level.dx = ", parent_level.dx)
                println("parent_level.x = ", parent_level.x)
                println("level_domain = ", level_domain)

                error(
                    "Level $(l): Coordinates are not aligned: ",
                    parent_level.x[parent_indices[1]],
                    " != ",
                    level_domain[1],
                    " or ",
                    parent_level.x[parent_indices[end]],
                    " != ",
                    level_domain[2],
                )
            end

            # build level
            level = Level{NumState,NumDiagnostic}(
                level_num_interior_points,
                num_ghost_points,
                num_buffer_points,
                num_transition_points,
                time_interpolation_order,
                spatial_interpolation_order,
                level_domain,
                physical_domain_box,
                level_dt,
                initial_time,
                false,
                parent_indices,
            )
            # ensure level is properly embedded in parent level
            if !(level.is_physical_boundary[1] || level.x[1] > parent_level.domain_box[1])
                error(
                    "Level $(l) is not properly embedded in parent level $(l-1), ",
                    "level.x[1] = $(level.x[1]), parent_level.domain_box[1] = $(parent_level.domain_box[1])",
                )
            end
            if !(level.is_physical_boundary[2] || level.x[end] < parent_level.domain_box[2])
                error(
                    "Level $(l) is not properly embedded in parent level $(l-1), ",
                    "level.x[end] = $(level.x[end]), parent_level.domain_box[2] = $(parent_level.domain_box[2])",
                )
            end
            push!(levels, level)
        end

        # construct
        return new{NumState,NumDiagnostic}(
            num_levels, levels, base_dt, initial_time, subcycling
        )
    end
end

"""
    Base.show(io::IO, grid::Grid)

Display a summary of the `Grid` structure, printing details for each level.
"""
function Base.show(
    io::IO, grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    println(io, "Grid Structure:")
    println(io, "  subcycling = ", grid.subcycling)
    for i in 1:(grid.num_levels)
        println(io, "level[", i, "],")
        println(io, "  num_interior_points       = ", grid.levels[i].num_interior_points)
        println(io, "  num_ghost_points          = ", grid.levels[i].num_ghost_points)
        println(io, "  num_buffer_points         = ", grid.levels[i].num_buffer_points)
        println(io, "  num_transition_points     = ", grid.levels[i].num_transition_points)
        println(
            io,
            "  spatial_interpolation_order = ",
            grid.levels[i].spatial_interpolation_order,
        )
        if length(grid.levels[i].parent_map) == grid.levels[i].num_total_points
            println(
                io,
                "  ibox   = ",
                [
                    grid.levels[i].parent_map[1 + grid.levels[i].num_buffer_points],
                    grid.levels[i].parent_map[grid.levels[i].num_total_points - grid.levels[i].num_buffer_points],
                ],
                ", ",
                [
                    grid.levels[i].is_aligned[1 + grid.levels[i].num_buffer_points],
                    grid.levels[i].is_aligned[grid.levels[i].num_total_points - grid.levels[i].num_buffer_points],
                ],
            )
        end
        println(io, "  domain_box  = ", grid.levels[i].domain_box)
        println(io, "  dx          = ", grid.levels[i].dx)
        println(io, "  dt          = ", grid.levels[i].dt)
        println(io, "  t           = ", grid.levels[i].t)
    end
end

"""
    shift_physical_boundary!(grid::Grid, num_shift_points::NTuple{2,Int})

Shift the physical boundary of the grid.

!!! warning
    This function is currently a placeholder and is not implemented.
"""
function shift_physical_boundary!(grid::Grid, num_shift_points::NTuple{2,Int}) end
