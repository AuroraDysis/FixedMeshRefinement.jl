export restriction!

#===============================================================================
restriction:
    * from level l+1 to level l
    * we assume that we always march fine level first (for l in lmax-1:-1:1)
    * we assume all the levels are at the same time slice
===============================================================================#
function restriction!(grid::Grid{NumState,NumDiagnostic}, l::Int; apply_trans_zone=false) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l + 1]
    coarse_level = grid.levels[l]

    (; num_total_points, num_left_ghost_points, num_right_ghost_points, num_transition_points, parent_indices) =
        fine_level

    isrt = if apply_trans_zone
        1 + num_buffer_points + num_transition_points
    else
        1 + num_buffer_points
    end

    iend = if apply_trans_zone
        num_total_points - num_buffer_points - num_transition_points
    else
        num_total_points - num_buffer_points
    end

    uf = view(fine_level.state[end], isrt:2:iend, :)
    uc = coarse_level.state[end]

    # TODO: implement fourth-order interpolation
    if apply_trans_zone
        noffset = div(num_transition_points + 1, 2)
        uc[parent_indices[(1 + noffset):(end - noffset)], :] .= uf
    else
        uc[parent_indices, :] .= uf
    end
end
