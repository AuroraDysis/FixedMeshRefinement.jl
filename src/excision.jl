export shift_level_boundaries!, shift_grid_boundaries!

const DEFAULT_FILL_EXTENDED_GRID_EXTRAPOLATION_ORDER = 4

function _extrapolate!(out, input, order::Int)
    length(input) == order || error("Length of input must be equal to order.")

    if order == 1
        @.. out = input[1]
    elseif order == 2
        @.. out = 2 * input[1] - input[2]
    elseif order == 3
        @.. out = 3 * input[1] - 3 * input[2] + input[3]
    elseif order == 4
        @.. out = 4 * input[1] - 6 * input[2] + 4 * input[3] - input[4]
    elseif order == 5
        @.. out = 5 * input[1] - 10 * input[2] + 10 * input[3] - 5 * input[4] + input[5]
    elseif order == 6
        @.. out =
            6 * input[1] - 15 * input[2] + 20 * input[3] - 15 * input[4] + 6 * input[5] -
            input[6]
    else
        error("Extrapolation order $order is not supported.")
    end
end

function fill_extended_grid_extrapolate!(
    state, extended_indices, direction::Symbol, order::Int
)
    if order < 0
        error("Extrapolation order must be non-negative.")
    end

    if isempty(extended_indices)
        return nothing
    end

    # This implementation assumes that `state` is an iterable over the arrays to be filled.
    # For example, a tuple of arrays: `(u, v, w)`. If `state` is a single array,
    # it should be wrapped in a tuple, e.g., `(state,)`.
    if direction == :left
        for idx in extended_indices
            input = [@view(state[i, :]) for i in (idx + 1):(idx + order)]
            _extrapolate!(@view(state[idx, :]), input, order)
        end
    elseif direction == :right
        for idx in extended_indices
            input = [@view(state[i, :]) for i in (idx - 1):-1:(idx - order)]
            _extrapolate!(@view(state[idx, :]), input, order)
        end
    else
        error("Unsupported direction for extrapolation: $direction")
    end

    return nothing
end

"""
    shift_level_boundaries!(level::Level, num_shift_points::NTuple{2,Int}; func_fill_extended::Function=(i, level) -> 0.0)

Shift the boundary of a `Level` by a given number of points, as a prerequisite, the level must align with the physical boundary.
"""
function shift_level_boundaries!(
    level::Level,
    num_shift_points::NTuple{2,Int};
    func_fill_extended::Function=(state, extended_indices, direction) ->
        fill_extended_grid_extrapolate!(
            state,
            extended_indices,
            direction,
            DEFAULT_FILL_EXTENDED_GRID_EXTRAPOLATION_ORDER,
        ),
)
    (; is_physical_boundary) = level

    if (num_shift_points[1] != 0 && !is_physical_boundary[1]) ||
        (num_shift_points[2] != 0 && !is_physical_boundary[2])
        println("level.domain_box = ", level.domain_box)
        println("level.is_physical_boundary = ", level.is_physical_boundary)
        error(
            "shifting the boundary of the level that is not aligned with the physical boundary is not allowed.",
        )
    end

    (;
        num_interior_points,
        domain_box,
        dx,
        physical_domain_box,
        parent_indices,
        offset_indices,
        is_base_level,
    ) = level

    # if base level then num_shift_points can be odd, otherwise it must be even
    if !is_base_level
        num_shift_points[1] % 2 == 0 && num_shift_points[2] % 2 == 0 || error(
            "num_shift_points must be even for non-base levels, num_shift_points = $(num_shift_points)",
        )
    end

    new_num_interior_points =
        num_interior_points + num_shift_points[1] + num_shift_points[2]
    new_domain_box = (
        domain_box[1] - num_shift_points[1] * dx, domain_box[2] + num_shift_points[2] * dx
    )
    new_physical_domain_box = (
        physical_domain_box[1] - num_shift_points[1] * dx,
        physical_domain_box[2] + num_shift_points[2] * dx,
    )
    new_parent_indices = if is_base_level
        0:0
    else
        left_parent_indices = first(parent_indices) - div(num_shift_points[1], 2)
        right_parent_indices = last(parent_indices) + div(num_shift_points[2], 2)
        left_parent_indices:right_parent_indices
    end
    new_offset_indices = offset_indices .+ num_shift_points[1]

    level.num_interior_points = new_num_interior_points
    level.domain_box = new_domain_box
    level.physical_domain_box = new_physical_domain_box
    level.parent_indices = new_parent_indices
    level.offset_indices = new_offset_indices

    # Fill extended regions
    if num_shift_points[1] > 0
        state = get_state(level)
        extended_indices = num_shift_points[1]:-1:1
        func_fill_extended(state, extended_indices, :left)
    end

    if num_shift_points[2] > 0
        state = get_state(level)
        extended_indices =
            (new_num_interior_points - num_shift_points[2] + 1):new_num_interior_points
        func_fill_extended(state, extended_indices, :right)
    end

    return nothing
end

"""
    shift_grid_boundaries!(grid::Grid, num_shift_points::NTuple{2,Int}; func_fill_extended::Function=(i, level) -> 0.0)

Shift the boundary of the grid by a given number of points; as a prerequisite, all levels must align with the physical boundary.
If the number of shift points is positive, the grid will be extended, otherwise it will be shrunk.
The tuple `num_shift_points` is the number of points to shift on the left and right boundaries, respectively.

# Arguments
- `grid::Grid`: The grid to shift boundaries of
- `num_shift_points::NTuple{2,Int}`: Number of points to shift on left and right boundaries
- `func_fill_extended::Function=(i, level) -> 0.0`: Function to compute fill values, takes index and level as arguments
"""
function shift_grid_boundaries!(
    grid::Grid,
    num_shift_points::NTuple{2,Int};
    func_fill_extended::Function=(state, extended_indices, direction) ->
        fill_extended_grid_extrapolate!(
            state,
            extended_indices,
            direction,
            DEFAULT_FILL_EXTENDED_GRID_EXTRAPOLATION_ORDER,
        ),
)
    (; num_levels, levels) = grid

    # make sure all levels are at same time
    t = levels[1].t
    for l in 2:num_levels
        level = levels[l]
        isapprox_tol(level.t, t) ||
            error("Level $(l) is not at same time, t = $(t), levels[l].t = $(levels[l].t)")
        if (num_shift_points[1] != 0 && !level.is_physical_boundary[1]) ||
            (num_shift_points[2] != 0 && !level.is_physical_boundary[2])
            error(
                "shifting the boundary of the level that is not aligned with the physical boundary is not allowed.",
            )
        end
    end

    for l in 1:num_levels
        level = levels[l]
        left_shift_points = num_shift_points[1] * 2^(l - 1)
        right_shift_points = num_shift_points[2] * 2^(l - 1)
        shift_level_boundaries!(
            level,
            (left_shift_points, right_shift_points);
            func_fill_extended=func_fill_extended,
        )
    end
end
