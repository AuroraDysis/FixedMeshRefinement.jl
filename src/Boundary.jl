module Boundary

function ApplyPeriodicBoundaryCondition!(
    grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    base_level = grid.levels[1]
    num_total_points = base_level.num_total_points
    num_buffer_points = base_level.num_buffer_points
    for v in 1:NumState
        u = base_level.state[v]
        for i in 1:num_buffer_points
            u[i] = u[num_total_points - 2 * num_buffer_points + i]
        end
        for i in num_total_points:-1:(num_total_points - num_buffer_points + 1)
            u[i] = u[2 * num_buffer_points - num_total_points + i]
        end
    end
end

function ApplyPeriodicBoundaryConditionRHS!(level, r)
    num_total_points = level.num_total_points
    num_buffer_points = level.num_buffer_points
    for v in 1:length(r)
        rhs = r[v]
        for i in 1:num_buffer_points
            rhs[i] = rhs[num_total_points - 2 * num_buffer_points + i]
        end
        for i in num_total_points:-1:(num_total_points - num_buffer_points + 1)
            rhs[i] = rhs[2 * num_buffer_points - num_total_points + i]
        end
    end
end

end
