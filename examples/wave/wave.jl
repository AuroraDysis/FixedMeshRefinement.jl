
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
            ) / dx
        diss_Pi =
            (
                (Pi[i + 3] + Pi[i - 3]) - 6 * (Pi[i + 2] + Pi[i - 2]) +
                15 * (Pi[i + 1] + Pi[i - 1]) - 20 * Pi[i]
            ) / dx
        psi_rhs[i] = Pi[i] + dissipation * diss_psi
        Pi_rhs[i] = ddpsi + dissipation * diss_Pi
    end

    psi_rhs[begin:(first(rhs_indices) - 1)] .= NaN
    psi_rhs[(last(rhs_indices) + 1):end] .= NaN

    apply_reflective_boundary_condition_rhs!(level, rhs)

    return nothing
end

#===============================================================================
Energy:
    * int_xmin^xmax (Pi^2/2 + dpsi^2/2)
    * calculate on base level (interior) only
===============================================================================#
function integrate(y::AbstractVector{Float64}, dx::Float64)
    @inbounds retval =
        (
            17 * (y[1] + y[end]) +
            59 * (y[2] + y[end - 1]) +
            43 * (y[3] + y[end - 2]) +
            49 * (y[4] + y[end - 3])
        ) / 48
    @simd for i in 5:(length(y) - 4)
        @inbounds retval += y[i]
    end
    return retval * dx
end

function wave_energy(grid)
    base_level = grid.levels[1]
    (; dx) = base_level

    u = get_state(base_level)
    psi = @view(u[:, 1])
    Pi = @view(u[:, 2])

    diag_state = get_diagnostic_state(base_level)
    rho = @view(diag_state[:, 1])

    interior_indices = get_interior_indices(base_level)

    for i in interior_indices
        # 4th order finite difference
        dpsi = (psi[i - 2] - 8 * psi[i - 1] + 8 * psi[i + 1] - psi[i + 2]) / (12 * dx)
        rho[i] = (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi * dpsi)
    end

    # integrate over the domain
    E = integrate(@view(rho[interior_indices]), dx)

    return E
end
