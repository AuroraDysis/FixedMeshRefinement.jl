export restrict_injection!

"""
    restrict_injection!(grid::Grid, l::Int; apply_trans_zone=false)

Perform restriction from a fine grid (`l+1`) to a coarse grid (`l`). This operation
transfers data from a finer grid to a coarser grid. Currently, it uses injection
by taking the values from the fine grid at the locations corresponding to the
coarse grid points.

# Arguments
- `grid::Grid`: The grid structure.
- `l::Int`: The coarse level index.
- `apply_trans_zone::Bool`: If `true`, the restriction is only applied outside of the
  transition zones. Defaults to `false`.
"""
function restrict_injection!(grid::Grid, l::Int; apply_trans_zone=false)
    fine_level = grid.levels[l + 1]
    coarse_level = grid.levels[l]

    (; num_interior_points, num_transition_points, parent_indices) = fine_level

    if apply_trans_zone
        fidx_transition_left_offset = mod(num_transition_points[1], 2) == 0 ? 0 : 1
        fidx_transition_right_offset = mod(num_transition_points[2], 2) == 0 ? 0 : 1
        fidx_start = 1 + num_transition_points[1] + fidx_transition_left_offset
        fidx_end =
            num_interior_points - num_transition_points[2] - fidx_transition_right_offset
    else
        fidx_start = 1
        fidx_end = num_interior_points
    end

    uf = view(get_state(fine_level), fidx_start:2:fidx_end, :)
    uc = get_state(coarse_level)

    if apply_trans_zone
        noffset_left = div(num_transition_points[1] + 1, 2)
        noffset_right = div(num_transition_points[2] + 1, 2)
        uc[parent_indices[(1 + noffset_left):(end - noffset_right)], :] .= uf
    else
        uc[parent_indices, :] .= uf
    end
end
