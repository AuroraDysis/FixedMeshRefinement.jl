function apply_reflective_boundary_condition!(
    grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    (; num_levels, levels) = grid

    for l in 1:num_levels
        level = levels[l]
        u = level.state[end]
        (; num_total_points, num_additional_points, additional_points_indices, is_physical_boundary) = level

        # apply reflective boundary condition to state
        if is_physical_boundary[1]
            for i in 1:num_additional_points[1]
                @. u[additional_points_indices[1][i], :] = -@view(u[num_additional_points[1] + i, :])
            end
        end

        if is_physical_boundary[2]
            for i in 1:num_additional_points[2]
                @. u[additional_points_indices[2][i], :] = -@view(
                    u[num_total_points - num_additional_points[2] + 1 - i, :]
                )
            end
        end
    end

    return nothing
end

function apply_reflective_boundary_condition_rhs!(
    level::Level{NumState,NumDiagnostic}, rhs
) where {NumState,NumDiagnostic}
    (; num_total_points, num_additional_points, additional_points_indices, is_physical_boundary) = level

    # apply reflective boundary condition to state
    if is_physical_boundary[1]
        for i in 1:num_additional_points[1]
            @. rhs[additional_points_indices[1][i], :] = -@view(rhs[num_additional_points[1] + i, :])
        end
    end

    if is_physical_boundary[2]
        for i in 1:num_additional_points[2]
            @. rhs[additional_points_indices[2][i], :] = -@view(
                rhs[num_total_points - num_additional_points[2] + 1 - i, :]
            )
        end
    end

    return nothing
end
