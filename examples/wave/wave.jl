#===============================================================================
wave_rhs!:
    * rhs of wave equation
        dot(psi) = Pi
        dot(Pi)  = ddpsi
===============================================================================#
function wave_rhs!(level, rhs, u, p, t)
    psi = @view(u[:, 1])
    Pi = @view(u[:, 2])
    psi_rhs = @view(rhs[:, 1])
    Pi_rhs = @view(rhs[:, 2])

    (; dx) = level
    (; dissipation) = p

    rhs_indices = get_rhs_evaluation_indices(level)

    @inbounds for i in rhs_indices
        # 4th order finite difference
        ddpsi =
            (-psi[i - 2] + 16 * psi[i - 1] - 30 * psi[i] + 16 * psi[i + 1] - psi[i + 2]) /
            (12 * dx^2)
        diss_psi =
            (
                (psi[i + 3] + psi[i - 3]) - 6 * (psi[i + 2] + psi[i - 2]) +
                15 * (psi[i + 1] + psi[i - 1]) - 20 * psi[i]
            ) / (64 * dx)
        diss_Pi =
            (
                (Pi[i + 3] + Pi[i - 3]) - 6 * (Pi[i + 2] + Pi[i - 2]) +
                15 * (Pi[i + 1] + Pi[i - 1]) - 20 * Pi[i]
            ) / (64 * dx)
        psi_rhs[i] = Pi[i] + dissipation * diss_psi
        Pi_rhs[i] = ddpsi + dissipation * diss_Pi
    end

    first_idx = first(rhs_indices)
    last_idx = last(rhs_indices)
    psi_rhs[begin:(first_idx - 1)] .= NaN
    psi_rhs[(last_idx + 1):end] .= NaN
    Pi_rhs[begin:(first_idx - 1)] .= NaN
    Pi_rhs[(last_idx + 1):end] .= NaN

    apply_reflective_boundary_condition_rhs!(level, rhs)

    return nothing
end

#===============================================================================
Energy:
    * int_xmin^xmax (Pi^2/2 + dpsi^2/2)
    * calculate on base level (interior) only
===============================================================================#
function wave_energy(grid)
    apply_reflective_boundary_condition!(grid)

    base_level = grid.levels[1]
    (; dx) = base_level

    u = get_state(base_level)
    psi = @view(u[:, 1])
    Pi = @view(u[:, 2])

    diag_state = get_diagnostic_state(base_level)
    rho = @view(diag_state[:, 1])

    interior_indices = get_interior_indices(base_level)

    @inbounds for i in interior_indices
        # 4th order finite difference
        dpsi = (psi[i - 2] - 8 * psi[i - 1] + 8 * psi[i + 1] - psi[i + 2]) / (12 * dx)
        rho[i] = (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi * dpsi)
    end

    # TODO: improve trapezoidal rule for boundary points
    first_idx = first(interior_indices)
    last_idx = last(interior_indices)
    indices = (first_idx + 1):(last_idx - 1)
    E =
        sum(@view(rho[indices])) * dx +
        (rho[first_idx] + rho[last_idx]) * dx / 2 +
        0.5 * (rho[first_idx] + rho[last_idx]) * dx

    return E
end
