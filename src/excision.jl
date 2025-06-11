export shift_level_boundaries!, shift_grid_boundaries!

"""
    shift_level_boundaries!(level::Level, num_shift_points::NTuple{2,Int})

Shift the boundary of a `Level` by a given number of points, as a prerequisite, the level must align with the physical boundary.
"""
function shift_level_boundaries!(level::Level, num_shift_points::NTuple{2,Int})
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
        num_interior_points - num_shift_points[1] - num_shift_points[2]
    new_domain_box = (
        domain_box[1] + num_shift_points[1] * dx, domain_box[2] - num_shift_points[2] * dx
    )
    new_physical_domain_box = (
        physical_domain_box[1] + num_shift_points[1] * dx,
        physical_domain_box[2] - num_shift_points[2] * dx,
    )
    new_parent_indices = if is_base_level
        0:0
    else
        left_parent_indices = first(parent_indices) + div(num_shift_points[1], 2)
        right_parent_indices = last(parent_indices) - div(num_shift_points[2], 2)
        left_parent_indices:right_parent_indices
    end
    new_offset_indices = offset_indices .- num_shift_points[1]

    level.num_interior_points = new_num_interior_points
    level.domain_box = new_domain_box
    level.physical_domain_box = new_physical_domain_box
    level.parent_indices = new_parent_indices
    level.offset_indices = new_offset_indices

    return nothing
end

"""
    shift_grid_boundaries!(grid::Grid, num_shift_points::NTuple{2,Int})

Shift the boundary of the grid by a given number of points; as a prerequisite, all levels must align with the physical boundary.
If the number of shift points is negative, the grid will be extended, otherwise it will be shrunk.
"""
function shift_grid_boundaries!(grid::Grid, num_shift_points::NTuple{2,Int})
    (; num_levels, levels) = grid

    # make sure all levels are at same time
    t = levels[1].t
    for l in 2:num_levels
        level = levels[l]
        isapprox_tol(level.t, t) ||
            error("Level $(l) is not at same time, t = $(t), levels[l].t = $(levels[l].t)")
    end

    for l in 1:num_levels
        level = levels[l]
        left_shift_points = num_shift_points[1] * 2^(l - 1)
        right_shift_points = num_shift_points[2] * 2^(l - 1)
        shift_level_boundaries!(level, (left_shift_points, right_shift_points))
    end
end
