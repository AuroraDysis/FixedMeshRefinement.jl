export excise_level!, excise_grid!

"""
    excise_level!(level::Level, num_excise_points::NTuple{2,Float64})

Excise the boundary of a `Level` by a given number of points, as a prerequisite, the level must align with the physical boundary.
"""
function excise_level!(level::Level, num_excise_points::NTuple{2,Float64})
    (; is_physical_boundary) = level

    if (num_excise_points[1] > 0 && !is_physical_boundary[1]) ||
        (num_excise_points[2] > 0 && !is_physical_boundary[2])
        println("level.domain_box = ", level.domain_box)
        println("level.is_physical_boundary = ", level.is_physical_boundary)
        error("Level is not aligned with the physical boundary")
    end

    (;
        num_interior_points,
        num_additional_points,
        num_ghost_points,
        domain_box,
        offset_indices,
        dx,
        physical_domain_box,
        parent_indices,
        is_base_level,
        state,
        rhs,
        tmp,
        k,
        diag_state,
    ) = level

    # if base level then num_excise_points can be odd, otherwise it must be even
    if !is_base_level
        num_excise_points[1] % 2 == 0 && num_excise_points[2] % 2 == 0 ||
            error("num_excise_points must be even for non-base levels")
    end

    new_num_interior_points =
        num_interior_points - num_excise_points[1] - num_excise_points[2]
    new_num_total_points =
        new_num_interior_points + num_additional_points[1] + num_additional_points[2]
    new_domain_box = (
        domain_box[1] + num_excise_points[1] * dx, domain_box[2] - num_excise_points[2] * dx
    )
    new_physical_domain_box = (
        physical_domain_box[1] + num_excise_points[1] * dx,
        physical_domain_box[2] - num_excise_points[2] * dx,
    )
    new_parent_indices = if is_base_level
        (0, 0)
    else
        parent_indices[(1 + div(num_excise_points[1], 2)):(end - div(num_excise_points[2], 2))]
    end
    new_additional_points_indices = (
        num_additional_points[1], # don't change the left indices since we use offset array
        num_additional_points[2] .- num_excise_points[2],
    )
    x_min = new_domain_box[1] - num_additional_points[1] * dx
    x_max = new_domain_box[2] + num_additional_points[2] * dx
    new_offset_indices = offset_indices .- num_excise_points[1]
    new_rhs_indices =
        (-num_additional_points[1] + 1 + num_ghost_points):(new_num_interior_points + num_additional_points[2] - num_ghost_points)
    new_x = OffsetVector(LinRange(x_min, x_max, new_num_total_points), new_offset_indices)
    new_state = [OffsetArray(s, new_offset_indices, :) for s in state]
    new_rhs = OffsetArray(rhs, new_offset_indices, :)
    new_tmp = OffsetArray(tmp, new_offset_indices, :)
    new_k = [OffsetArray(ki, new_offset_indices, :) for ki in k]
    new_diag_state = OffsetArray(diag_state, new_offset_indices, :)

    # TODO: assign new values to level
    level.num_interior_points = new_num_interior_points
    level.num_total_points = new_num_total_points
    level.domain_box = new_domain_box
    level.physical_domain_box = new_physical_domain_box
    level.parent_indices = new_parent_indices
    level.additional_points_indices = new_additional_points_indices
    level.rhs_indices = new_rhs_indices
    level.x = new_x
    level.state = new_state
    level.rhs = new_rhs
    level.tmp = new_tmp
    level.k = new_k
    level.diag_state = new_diag_state

    return nothing
end

"""
    excise_grid!(grid::Grid, num_excise_points::NTuple{2,Float64})

Excise the grid by a given number of points; as a prerequisite, all levels must align with the physical boundary.
"""
function excise_grid!(grid::Grid, num_excise_points::NTuple{2,Float64})
    num_excise_points[1] >= 0 && num_excise_points[2] >= 0 ||
        error("num_excise_points must be non-negative")

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
        left_excise_points = num_excise_points[1] * 2^(l - 1)
        right_excise_points = num_excise_points[2] * 2^(l - 1)
        excise_level!(level, (left_excise_points, right_excise_points))
    end
end
