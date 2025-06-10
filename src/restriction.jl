export restriction!

"""
    restriction!(grid::Grid, l::Int; apply_trans_zone=false)

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
function restriction!(
    grid::Grid{NumState,NumDiagnostic}, l::Int; apply_trans_zone=false
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l + 1]
    coarse_level = grid.levels[l]

    (; num_interior_points, num_transition_points, parent_indices) = fine_level

    isrt = if apply_trans_zone
        1 + num_transition_points
    else
        1
    end

    iend = if apply_trans_zone
        num_interior_points - num_transition_points
    else
        num_interior_points
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
