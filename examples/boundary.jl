function apply_periodic_boundary_condition!(
    grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    base_level = grid.levels[1]
    (; num_total_points, num_buffer_points) = base_level

    u = base_level.state[end]

    for i in 1:num_buffer_points
        u[i, :] .= @view(u[num_total_points - 2 * num_buffer_points + i, :])
    end

    for i in (num_total_points - num_buffer_points + 1):num_total_points
        u[i, :] .= @view(u[2 * num_buffer_points - num_total_points + i, :])
    end
end

function apply_periodic_boundary_condition_rhs!(level, rhs)
    (; num_total_points, num_buffer_points) = level

    for i in 1:num_buffer_points
        rhs[i, :] .= @view(rhs[num_total_points - 2 * num_buffer_points + i, :])
    end

    for i in (num_total_points - num_buffer_points + 1):num_total_points
        rhs[i, :] .= @view(rhs[2 * num_buffer_points - num_total_points + i, :])
    end
end
