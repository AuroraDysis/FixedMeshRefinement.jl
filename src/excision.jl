export shift_level_boundaries!, shift_grid_boundaries!

const DEFAULT_FILL_EXTENDED_GRID_EXTRAPOLATION_ORDER = 5

function extrapolate!(out, input, order::Int)
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

    if direction == :left
        for idx in extended_indices
            input = [@view(state[i, :]) for i in (idx + 1):(idx + order)]
            extrapolate!(@view(state[idx, :]), input, order)
        end
    elseif direction == :right
        for idx in extended_indices
            input = [@view(state[i, :]) for i in (idx - 1):-1:(idx - order)]
            extrapolate!(@view(state[idx, :]), input, order)
        end
    else
        error("Unsupported direction for extrapolation: $direction")
    end

    return nothing
end

"""
    shift_level_boundaries!(level::Level, num_shift_points::NTuple{2,Int}; func_fill_extended::Function=fill_extended_grid_extrapolate!)

Shift the boundary of a `Level` by a given number of points, as a prerequisite, the level must align with the physical boundary.
"""
function shift_level_boundaries!(
    grid::Grid,
    l::Int,
    num_shift_points::NTuple{2,Int};
    shift_parent_indices::Bool=false,
    func_fill_extended::Function=(state, extended_indices, direction) ->
        fill_extended_grid_extrapolate!(
            state,
            extended_indices,
            direction,
            DEFAULT_FILL_EXTENDED_GRID_EXTRAPOLATION_ORDER,
        ),
)
    level = get_level(grid, l)

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
    ) = level

    # if base level then num_shift_points can be odd, otherwise it must be even
    if !is_base_level(level)
        num_shift_points[1] % 2 == 0 && num_shift_points[2] % 2 == 0 || error(
            "num_shift_points must be even for non-base levels, num_shift_points = $(num_shift_points)",
        )
    end

    new_num_unused_points = level.num_unused_points .- num_shift_points
    if any(new_num_unused_points .< 0)
        error(
            "num_unused_points must be non-negative, please consider increasing preallocation, num_unused_points = $(level.num_unused_points), num_shift_points = $(num_shift_points)",
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
    new_parent_indices = if is_base_level(level)
        0:0
    elseif shift_parent_indices
        left_parent_indices = first(parent_indices)
        right_parent_indices =
            last(parent_indices) + div(num_shift_points[2], 2) + div(num_shift_points[1], 2)
        left_parent_indices:right_parent_indices
    else
        left_parent_indices = first(parent_indices) - div(num_shift_points[1], 2)
        right_parent_indices = last(parent_indices) + div(num_shift_points[2], 2)
        left_parent_indices:right_parent_indices
    end
    new_offset_indices = offset_indices .+ num_shift_points[1]

    level.num_unused_points = new_num_unused_points
    level.num_interior_points = new_num_interior_points
    level.domain_box = new_domain_box
    level.physical_domain_box = new_physical_domain_box
    level.parent_indices = new_parent_indices
    level.offset_indices = new_offset_indices

    # check if the parent indices are valid
    if l > 1
        parent_level = get_level(grid, l - 1)
        parent_x = get_x(parent_level)
        level_domain = level.domain_box

        if !(
            isapprox_tol(parent_x[new_parent_indices[1]], level_domain[1]) &&
            isapprox_tol(parent_x[new_parent_indices[end]], level_domain[2])
        )
            println("parent_indices = ", new_parent_indices)
            println("parent_level.dx = ", parent_level.dx)
            println("parent_level.x = ", parent_x)
            println("level_domain = ", level_domain)

            error(
                "Level $(l): Coordinates are not aligned: ",
                parent_x[new_parent_indices[1]],
                " != ",
                level_domain[1],
                " or ",
                parent_x[new_parent_indices[end]],
                " != ",
                level_domain[2],
            )
        end
    end

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
    shift_grid_boundaries!(grid::Grid, num_shift_points::NTuple{2,Int}; func_fill_extended::Function=fill_extended_grid_extrapolate!)

Shift the boundary of the grid by a given number of points; as a prerequisite, all levels must align with the physical boundary.
If the number of shift points is positive, the grid will be extended, otherwise it will be shrunk.
The tuple `num_shift_points` is the number of points to shift on the left and right boundaries, respectively.

# Arguments
- `grid::Grid`: The grid to shift boundaries of
- `num_shift_points::NTuple{2,Int}`: Number of points to shift on left and right boundaries
- `func_fill_extended::Function=fill_extended_grid_extrapolate!`: Function to compute fill values, takes index and level as arguments
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
    num_levels = get_num_levels(grid)

    for l in 2:num_levels
        level = get_level(grid, l)
        if (num_shift_points[1] != 0 && !level.is_physical_boundary[1]) ||
            (num_shift_points[2] != 0 && !level.is_physical_boundary[2])
            error(
                "shifting the boundary of the level that is not aligned with the physical boundary is not allowed.",
            )
        end
    end

    for l in 1:num_levels
        level = get_level(grid, l)
        left_shift_points = num_shift_points[1] * 2^(l - 1)
        right_shift_points = num_shift_points[2] * 2^(l - 1)
        shift_level_boundaries!(
            grid,
            l,
            (left_shift_points, right_shift_points);
            shift_parent_indices=true,
            func_fill_extended=func_fill_extended,
        )
    end
end
