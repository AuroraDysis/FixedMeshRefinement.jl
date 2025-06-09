
#===============================================================================
wave_rhs!:
    * rhs of wave equation
        dot(psi) = Pi
        dot(Pi)  = ddpsi
===============================================================================#
function wave_rhs!(level, rhs, u, t)
    psi = u[1]
    Pi = u[2]
    psi_rhs = rhs[1]
    Pi_rhs = rhs[2]

    (; num_total_points, finite_difference_order, dissipation) = level

    # TODO: improve performance by using pre-allocated arrays
    ddpsi = zeros(Float64, num_total_points)
    psi_diss = zeros(Float64, num_total_points)
    Pi_diss = zeros(Float64, num_total_points)
    derivs_2nd!(ddpsi, psi, dx, finite_difference_order)
    derivs_diss!(psi_diss, psi, dx, finite_difference_order)
    derivs_diss!(Pi_diss, Pi, dx, finite_difference_order)

    @. psi_rhs = Pi + dissipation * psi_diss
    @. Pi_rhs = ddpsi + dissipation * Pi_diss

    if level.is_base_level
        Boundary.ApplyPeriodicBoundaryConditionRHS!(level, rhs)
    end
end

#===============================================================================
Energy:
    * int_xmin^xmax (Pi^2/2 + dpsi^2/2)
    * calculate on base level (interior) only
===============================================================================#
function wave_energy(grid)
    level = grid.levels[1]
    num_total_points = level.num_total_points
    num_buffer_points = level.num_buffer_points
    dx = level.dx
    u = level.state[end]
    psi = @view(u[:, 1])
    psi_t = @view(u[:, 2])

    E::Float64 = 0.0
    for i in (1 + num_buffer_points):(num_total_points - num_buffer_points)
        # 4th order finite difference
        dpsi = (-psi[i - 2] + 8 * psi[i - 1] - 8 * psi[i + 1] + psi[i + 2]) / (12 * dx)
        E += (0.5 * psi_t[i] * psi_t[i] + 0.5 * dpsi[i] * dpsi[i])
    end
    return E * dx
end
