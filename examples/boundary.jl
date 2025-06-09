function apply_reflective_boundary_condition!(
    grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    (; num_levels, levels) = grid

    for l in 1:num_levels
        level = levels[l]
        u = level.state[end]

        # apply reflective boundary condition to state
        if level.is_physical_boundary[1]
            for i in 1:level.num_ghost_points[1]
                u[level.ghost_indices[1][i], :] .= @view(
                    u[level.num_ghost_points[1] + i, :]
                )
            end
        end

        if level.is_physical_boundary[2]
            for i in 1:level.num_ghost_points[2]
                u[level.ghost_indices[2][i], :] .= @view(
                    u[level.num_total_points - level.num_ghost_points[2] + 1 - i, :]
                )
            end
        end
    end
end

function apply_reflective_boundary_condition_rhs!(level, rhs)
    (; num_total_points, num_buffer_points) = level

    for i in 1:num_buffer_points
        rhs[i, :] .= @view(rhs[num_total_points - 2 * num_buffer_points + i, :])
    end

    for i in (num_total_points - num_buffer_points + 1):num_total_points
        rhs[i, :] .= @view(rhs[2 * num_buffer_points - num_total_points + i, :])
    end
end
