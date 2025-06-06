module Boundary

function ApplyPeriodicBoundaryCondition!(gfs)
    lev1fs = grid.levels[1]
    num_total_points = grid.levels[1].num_total_points
    num_buffer_points = grid.levels[1].num_buffer_points
    for v = 1:grid.NumState
        u = lev1fs.u[v]
        for i = 1:num_buffer_points
            u[i] = u[num_total_points-2*num_buffer_points+i]
        end
        for i = num_total_points:-1:num_total_points-num_buffer_points+1
            u[i] = u[2*num_buffer_points-num_total_points+i]
        end
    end
end

function ApplyPeriodicBoundaryConditionRHS!(level, r)
    num_total_points = level.num_total_points
    num_buffer_points = level.num_buffer_points
    for v = 1:length(r)
        rhs = r[v]
        for i = 1:num_buffer_points
            rhs[i] = rhs[num_total_points-2*num_buffer_points+i]
        end
        for i = num_total_points:-1:num_total_points-num_buffer_points+1
            rhs[i] = rhs[2*num_buffer_points-num_total_points+i]
        end
    end
end

end
