export Level,
    Grid,
    get_x,
    get_state,
    get_rk4_tmp_state,
    get_rk_stage,
    get_diagnostic_state,
    get_boundary_indices,
    get_interior_indices,
    get_offset_indices,
    get_rhs_evaluation_indices,
    get_total_grid_points,
    get_maximum_grid_points,
    cycle_state!,
    fine_to_coarse_index

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
- `index::Int`: Index of the level.
- `num_interior_points::Int`: Number of interior grid points.
- `num_ghost_points::Int`: Number of ghost points on each side.
- `num_buffer_points::Int`: Number of buffer points on each side for inter-level communication.
- `num_transition_points::NTuple{2,Int}`: Number of points in the transition zone for mesh refinement.
- `num_boundary_points::NTuple{2,Int}`: Number of boundary points on each side (ghost or buffer).
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
- `offset_indices::UnitRange{Int}`: Range of indices for `OffsetArray`s of this level.
"""
mutable struct Level{NumState,NumDiagnostic}
    const index::Int
    num_interior_points::Int  # num of interior grid points
    const num_ghost_points::Int  # num of ghost points on each side
    const num_buffer_points::Int # num of buffer points on each side
    const num_transition_points::NTuple{2,Int}  # num of transition zone points
    const num_boundary_points::NTuple{2,Int} # num of boundary points on each side
    num_unused_points::NTuple{2,Int}
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
    offset_indices::UnitRange{Int}

    # data
    const x::LinRange{Float64,Int64}
    const state::Vector{Matrix{Float64}} # state vectors at different time levels
    const rk4_tmp_state::Matrix{Float64} # temporary state for RK4

    # intermediate state vectors for new subcycling
    const rk_stages::Vector{Matrix{Float64}}
    Yn_buffer::Vector{Array{Float64,3}}

    # diagnostic variables
    const diagnostic_state::Matrix{Float64} # state vectors for diagnostic variables

    function Level{NumState,NumDiagnostic}(
        index,
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
            isapprox_tol(domain_box[1], physical_domain_box[1]),
            isapprox_tol(domain_box[2], physical_domain_box[2]),
        )
        num_left_boundary_points =
            is_physical_boundary[1] ? num_ghost_points : num_buffer_points
        num_right_boundary_points =
            is_physical_boundary[2] ? num_ghost_points : num_buffer_points
        num_total_points =
            num_interior_points + num_left_boundary_points + num_right_boundary_points
        dx = (domain_box[2] - domain_box[1]) / (num_interior_points - 1)
        x_min = domain_box[1] - num_left_boundary_points * dx
        x_max = domain_box[2] + num_right_boundary_points * dx
        x = LinRange(x_min, x_max, num_total_points)
        time_levels = max(time_interpolation_order + 1, 2)
        state = [fill(NaN, num_total_points, NumState) for _ in 1:time_levels]
        tmp = fill(NaN, num_total_points, NumState)
        k = Vector{Matrix{Float64}}(undef, 4)
        Yn_buffer = Vector{Array{Float64,3}}(undef, 4)
        diag_state = fill(NaN, num_total_points, NumDiagnostic)
        for j in 1:4
            k[j] = fill(NaN, num_total_points, NumState)
            Yn_buffer[j] = fill(NaN, num_buffer_points, NumState, 2)
        end
        offset_indices =
            (-num_left_boundary_points + 1):(num_interior_points + num_right_boundary_points)
        num_left_transition_points = is_physical_boundary[1] ? 0 : num_transition_points
        num_right_transition_points = is_physical_boundary[2] ? 0 : num_transition_points

        return new{NumState,NumDiagnostic}(
            index,
            num_interior_points,
            num_ghost_points,
            num_buffer_points,
            (num_left_transition_points, num_right_transition_points),
            (num_left_boundary_points, num_right_boundary_points),
            (0, 0),
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
            offset_indices,
            # data
            x,
            state,
            tmp,
            k,
            Yn_buffer,
            diag_state,
        )
    end
end

"""
    get_rhs_evaluation_indices(level::Level) -> UnitRange{Int}

Return the indices that requires evaluation of the right-hand side.
"""
function get_rhs_evaluation_indices(level::Level)
    (; num_boundary_points, num_interior_points, num_ghost_points) = level
    return (-num_boundary_points[1] + 1 + num_ghost_points):(num_interior_points + num_boundary_points[2] - num_ghost_points)
end

"""
    get_interior_indices(level::Level) -> UnitRange{Int}

Return the indices of the interior grid points.
"""
function get_interior_indices(level::Level)
    (; num_interior_points) = level
    return 1:num_interior_points
end

"""
    get_boundary_indices(level::Level) -> NTuple{2,StepRange{Int,Int}}

Return the indices of the boundary points on each side (ghost or buffer).
"""
function get_boundary_indices(level::Level)
    (; num_boundary_points, num_interior_points) = level
    return (
        0:-1:(-num_boundary_points[1] + 1),
        (num_interior_points + 1):(num_interior_points + num_boundary_points[2]),
    )
end

"""
    get_offset_indices(level::Level) -> UnitRange{Int}

Return the indices of the grid points in the `OffsetArray` of the `level`.
"""
function get_offset_indices(level::Level)
    return level.offset_indices
end

"""
    get_total_grid_points(level::Level) -> Int

Return the total number of grid points at this level (including ghost and buffer points).
"""
function get_total_grid_points(level::Level)
    (; num_interior_points, num_boundary_points) = level
    return num_interior_points + num_boundary_points[1] + num_boundary_points[2]
end

"""
    get_maximum_grid_points(level::Level) -> Int

Return the maximum number of grid points at this level (including ghost, buffer and unused points).
"""
function get_maximum_grid_points(level::Level)
    return length(level.x)
end

"""
    get_x(level::Level)

Return the grid point coordinates of the `level` as an `OffsetArray`.
"""
function get_x(level::Level)
    (; x, offset_indices) = level
    return OffsetArray(x, offset_indices)
end

"""
    get_state(level::Level, i::Int=0)

Return the state variables of the `level` as an `OffsetArray`. The optional
argument `i` specifies the time level, where `i=0` corresponds to the current
time level `n`, `i=-1` to `n-1`, etc.
"""
function get_state(level::Level, i::Int=0)
    (; state, offset_indices) = level
    return OffsetArray(state[end + i], offset_indices, :)
end

"""
    get_rk4_tmp_state(level::Level)

Return a temporary array with the same size as the state variables of the `level`
as an `OffsetArray`.
"""
function get_rk4_tmp_state(level::Level)
    (; rk4_tmp_state, offset_indices) = level
    return OffsetArray(rk4_tmp_state, offset_indices, :)
end

"""
    get_rk_stage(level::Level, i::Int)

Return the `i`-th intermediate state `k_i` for the Runge-Kutta time stepping
scheme as an `OffsetArray`.
"""
function get_rk_stage(level::Level, i::Int)
    (; rk_stages, offset_indices) = level
    return OffsetArray(rk_stages[i], offset_indices, :)
end

"""
    get_diagnostic_state(level::Level)

Return the diagnostic state variables of the `level` as an `OffsetArray`.
"""
function get_diagnostic_state(level::Level)
    (; diagnostic_state, offset_indices) = level
    return OffsetArray(diagnostic_state, offset_indices, :)
end

"""
    cycle_state!(level::Level)

Shift the time levels of the state variables in a `Level`. `state[i]` becomes
`state[i+1]`. This is used to advance the solution in time.
"""
function cycle_state!(level::Level)
    (; state) = level
    # shift state vectors
    for i in 1:(length(state) - 1)
        state[i] .= state[i + 1]
    end
    return nothing
end

"""
    fine_to_coarse_index(fine_level::Level, fidx::Int) -> Int

Convert a fine grid index `fidx` to the corresponding coarse grid index.
This function is not defined for the base level. An error is thrown if the fine
grid point does not align with any coarse grid point.
"""
function fine_to_coarse_index(fine_level::Level, fidx::Int)
    (; is_base_level, parent_indices) = fine_level

    if is_base_level
        error("fine_to_coarse_index is not defined for base level")
    end

    parent_idx_left = parent_indices[1]
    offset = fidx - 1

    if !iseven(offset)
        println("parent_indices = ", parent_indices)
        println("fidx = $fidx, offset = $offset")
        error("fidx = $fidx is not aligned with any point in parent level")
    end

    half_offset = div(offset, 2)
    cidx = parent_idx_left + half_offset
    return cidx
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

    """
        Grid{NumState,NumDiagnostic}(
            base_level_num_points,
            domain_boxes,
            num_ghost_points,
            num_buffer_points;
            kwargs...
        )

    Construct a `Grid`.

    # Arguments
    - `base_level_num_points::Int`: Number of interior grid points at the base level.
    - `domain_boxes::Vector{Tuple{Float64,Float64}}`: A vector of domain boxes for each level. `domain_boxes[1]` is the physical domain.
    - `num_ghost_points::Int`: Number of ghost points on each side.
    - `num_buffer_points::Int`: Number of buffer points for inter-level communication.

    # Keyword Arguments
    - `num_transition_points::Int=3`: Number of points in the transition zone for mesh refinement.
    - `time_interpolation_order::Int=2`: Order of time interpolation.
    - `spatial_interpolation_order::Int=5`: Order of spatial interpolation.
    - `cfl::Float64=0.25`: CFL number for time step calculation.
    - `initial_time::Float64=0.0`: Initial time of the simulation.
    - `subcycling::Bool=true`: Enable subcycling in time.
    """
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
            1,
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

            parent_x = get_x(parent_level)
            # check x are aligned
            if !(
                isapprox_tol(parent_x[parent_indices[1]], level_domain[1]) &&
                isapprox_tol(parent_x[parent_indices[end]], level_domain[2])
            )
                println("parent_indices = ", parent_indices)
                println("parent_level.dx = ", parent_level.dx)
                println("parent_level.x = ", parent_x)
                println("level_domain = ", level_domain)

                error(
                    "Level $(l): Coordinates are not aligned: ",
                    parent_x[parent_indices[1]],
                    " != ",
                    level_domain[1],
                    " or ",
                    parent_x[parent_indices[end]],
                    " != ",
                    level_domain[2],
                )
            end

            # build level
            level = Level{NumState,NumDiagnostic}(
                l,
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
            level_x = get_x(level)
            if !(level.is_physical_boundary[1] || level_x[1] > parent_level.domain_box[1])
                error(
                    "Level $(l) is not properly embedded in parent level $(l-1), ",
                    "level.x[1] = $(level_x[1]), parent_level.domain_box[1] = $(parent_level.domain_box[1])",
                )
            end
            if !(level.is_physical_boundary[2] || level_x[end] < parent_level.domain_box[2])
                error(
                    "Level $(l) is not properly embedded in parent level $(l-1), ",
                    "level.x[end] = $(level_x[end]), parent_level.domain_box[2] = $(parent_level.domain_box[2])",
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

Display a compact summary of the `Grid`.
"""
function Base.show(
    io::IO, grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    return print(
        io, "Grid{$NumState, $NumDiagnostic} with $(grid.num_levels) levels at t=$(grid.t)"
    )
end

"""
    Base.show(io::IO, ::MIME"text/plain", grid::Grid)

Display a detailed summary of the `Grid` structure.
"""
function Base.show(
    io::IO, ::MIME"text/plain", grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    println(
        io, "Grid{$NumState, $NumDiagnostic} with $(grid.num_levels) levels at t=$(grid.t):"
    )
    println(io, "  Subcycling: ", grid.subcycling)
    println(io, "  Base dt:    ", grid.base_dt)
    println(io, "Levels:")
    for i in 1:(grid.num_levels)
        println(io, "  [$i]: ", grid.levels[i])
    end
end

"""
    Base.show(io::IO, level::Level)

Display a compact summary of the `Level`.
"""
function Base.show(io::IO, level::Level)
    return print(
        io,
        "Level(index=",
        level.index,
        ", domain=",
        level.domain_box,
        ", N=",
        level.num_interior_points,
        ", dx=",
        round(level.dx; sigdigits=3),
        ")",
    )
end

"""
    Base.show(io::IO, ::MIME"text/plain", level::Level)

Display a detailed summary of the `Level`.
"""
function Base.show(
    io::IO, ::MIME"text/plain", level::Level{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    println(io, "Level{$NumState, $NumDiagnostic}:")
    println(io, "  Index:                 ", level.index)
    println(io, "  Domain:                ", level.domain_box)
    println(io, "  Interior points:       ", level.num_interior_points)
    println(io, "  Grid spacing (dx):     ", level.dx)
    println(io, "  Time step (dt):        ", level.dt)
    println(io, "  Current time (t):      ", level.t)
    println(io, "  Ghost points:          ", level.num_ghost_points)
    println(io, "  Buffer points:         ", level.num_buffer_points)
    println(io, "  Boundary points:       ", level.num_boundary_points)
    println(io, "  Transition points:     ", level.num_transition_points)
    println(io, "  Interpolation (time):  ", level.time_interpolation_order)
    println(io, "  Interpolation (space): ", level.spatial_interpolation_order)
    println(io, "  Physical boundary:     ", level.is_physical_boundary)
    println(io, "  Base level:            ", level.is_base_level)
    if !level.is_base_level
        println(io, "  Parent indices:        ", level.parent_indices)
    end
    return println(io, "  Indices:               ", level.offset_indices)
end
