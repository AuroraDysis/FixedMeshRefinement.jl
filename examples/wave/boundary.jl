function apply_reflective_boundary_condition!(
    grid::Grid{NumState,NumDiagnostic,NumTemp}
) where {NumState,NumDiagnostic,NumTemp}
    (; num_levels, levels) = grid

    for l in 1:num_levels
        level = levels[l]
        u = get_state(level)
        (; num_interior_points, num_boundary_points, is_physical_boundary) = level

        boundary_indices = get_boundary_indices(level)

        # apply reflective boundary condition to state
        if is_physical_boundary[1]
            for i in 1:num_boundary_points[1]
                @.. u[boundary_indices[1][i], :] = -@view(u[1 + i, :])
            end
        end

        if is_physical_boundary[2]
            for i in 1:num_boundary_points[2]
                @.. u[boundary_indices[2][i], :] = -@view(u[num_interior_points - i, :])
            end
        end
    end

    return nothing
end

function apply_reflective_boundary_condition_rhs!(
    level::Level{NumState,NumDiagnostic,NumTemp}, rhs
) where {NumState,NumDiagnostic,NumTemp}
    (; num_interior_points, num_boundary_points, is_physical_boundary) = level

    boundary_indices = get_boundary_indices(level)

    # apply reflective boundary condition to state
    if is_physical_boundary[1]
        for i in 1:num_boundary_points[1]
            @.. rhs[boundary_indices[1][i], :] = -@view(rhs[1 + i, :])
        end
    end

    if is_physical_boundary[2]
        for i in 1:num_boundary_points[2]
            @.. rhs[boundary_indices[2][i], :] = -@view(rhs[num_interior_points - i, :])
        end
    end

    return nothing
end
