
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

    (; num_ghost_points, x, dx) = level
    (; dissipation) = p

    @inbounds for i in (first(x) + num_ghost_points):(last(x) - num_ghost_points)
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

    psi_rhs[begin:0] .= NaN
    psi_rhs[(num_interior_points + 1):end] .= NaN

    return apply_reflective_boundary_condition_rhs!(level, rhs)
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
    (; num_interior_points, dx, state, diag_state) = base_level

    u = state[end]
    psi = @view(u[:, 1])
    Pi = @view(u[:, 2])

    rho = @view(diag_state[:, 1])

    for i in 1:num_interior_points
        # 4th order finite difference
        dpsi = (psi[i - 2] - 8 * psi[i - 1] + 8 * psi[i + 1] - psi[i + 2]) / (12 * dx)
        rho[i] = (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi * dpsi)
    end

    # integrate over the domain
    E = integrate(@view(rho[1:num_interior_points]), dx)

    return E
end
