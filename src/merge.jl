export merge_grid_levels

"""
    merge_grid_levels(grid::Grid, getter::Function)

Merge signle variable from all levels of the grid into a `BlockArray`.
The `getter` function specifies which data to retrieve the data of the variable from each level.
`getter(level)` should return an `OffsetArray` from which data is read.
"""
function merge_grid_levels(grid::Grid, getter::Function)
    num_levels = get_num_levels(grid)

    # Start with the interior of the finest level
    finest_level = get_level(grid, num_levels)
    finest_x = get_x(finest_level)
    finest_data = getter(finest_level)
    finest_interior_indices = get_interior_indices(finest_level)
    x_blocks = [@view finest_x[finest_interior_indices]]
    y_blocks = [@view finest_data[finest_interior_indices]]

    # Go from finest - 1 to coarsest, adding the uncovered parts of each level
    for l in (num_levels - 1):-1:1
        fine_level = get_level(grid, l + 1)
        level = get_level(grid, l)
        level_x = get_x(level)
        level_data = getter(level)
        level_interior_indices = get_interior_indices(level)
        l_lo_idx, l_hi_idx = first(level_interior_indices), last(level_interior_indices)
        f_lo_idx, f_hi_idx = first(fine_level.parent_indices),
        last(fine_level.parent_indices)

        # Prepend the part of the parent grid before the fine grid
        if l_lo_idx < f_lo_idx
            pushfirst!(x_blocks, @view level_x[l_lo_idx:(f_lo_idx - 1)])
            pushfirst!(y_blocks, @view level_data[l_lo_idx:(f_lo_idx - 1)])
        end

        # Append the part of the coarse grid after the fine grid
        if l_hi_idx > f_hi_idx
            push!(x_blocks, @view level_x[(f_hi_idx + 1):l_hi_idx])
            push!(y_blocks, @view level_data[(f_hi_idx + 1):l_hi_idx])
        end
    end

    return mortar(x_blocks), mortar(y_blocks)
end
