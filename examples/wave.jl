
#===============================================================================
wave_rhs!:
    * rhs of wave equation
        dot(psi) = Pi
        dot(Pi)  = ddpsi
===============================================================================#
function wave_rhs!(level, rhs, u, t)
    psi = @view(u[:, 1])
    Pi = @view(u[:, 2])
    psi_rhs = @view(rhs[:, 1])
    Pi_rhs = @view(rhs[:, 2])

    (;
        is_base_level,
        num_total_points,
        num_ghost_points,
        num_buffer_points,
        dissipation,
        dx,
    ) = level

    noffset = is_base_level ? num_ghost_points : num_buffer_points
    @inbounds for i in (1 + noffset):(num_total_points - noffset)
        # 4th order finite difference
        ddpsi =
            (-psi[i - 2] + 16 * psi[i - 1] - 30 * psi[i] + 16 * psi[i + 1] - psi[i + 2]) /
            (12 * dx^2)
        diss_psi =
            (
                (psi[i + 3] + psi[i - 3]) - 6 * (psi[i + 2] + psi[i - 2]) +
                15 * (psi[i + 1] + psi[i - 1]) - 20 * psi[i]
            ) / dx
        diss_Pi =
            (
                (Pi[i + 3] + Pi[i - 3]) - 6 * (Pi[i + 2] + Pi[i - 2]) +
                15 * (Pi[i + 1] + Pi[i - 1]) - 20 * Pi[i]
            ) / dx
        psi_rhs[i] = Pi[i] + dissipation * diss_psi
        Pi_rhs[i] = ddpsi + dissipation * diss_Pi
    end

    if is_base_level
        apply_periodic_boundary_condition_rhs!(level, rhs)
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
    Pi = @view(u[:, 2])

    E::Float64 = 0.0
    for i in (1 + num_buffer_points):(num_total_points - num_buffer_points)
        # 4th order finite difference
        dpsi = (-psi[i - 2] + 8 * psi[i - 1] - 8 * psi[i + 1] + psi[i + 2]) / (12 * dx)
        E += (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi * dpsi)
    end
    return E * dx
end
