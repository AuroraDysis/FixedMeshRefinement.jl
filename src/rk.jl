function rk4!(f::Function, level)
    (; u, u_p, rhs, tmp, t, dt) = level
    sixth_dt = dt / 6
    third_dt = dt / 3
    half_dt = dt / 2

    cycle_state!(level)

    level.t = t
    f(level, rhs, u)
    @. u += sixth_dt * rhs

    @. tmp = u_p + half_dt * rhs
    level.t = t + half_dt
    f(level, rhs, tmp)
    @. u += third_dt * rhs

    @. tmp = u_p + half_dt * rhs
    level.t = t + half_dt
    f(level, rhs, tmp)
    @. u += third_dt * rhs

    @. tmp = u_p + dt * rhs
    level.t = t + dt
    f(level, rhs, tmp)
    @. u += sixth_dt * rhs
    return level.t = t + dt
end

function rk4_mongwane!(
    f::Function, level::Level{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    (; u, u_p, rhs, tmp, k, t, dt) = level
    sixth_dt = dt / 6
    third_dt = dt / 3
    half_dt = dt / 2

    k1, k2, k3, k4 = k

    isrt = level.is_base_level ? 1 : 1 + level.num_buffer_points
    iend = if level.is_base_level
        level.num_total_points
    else
        level.num_total_points - level.num_buffer_points
    end

    @. u_p = u
    level.t = t
    f(level, rhs, u)
    k1[isrt:iend, :] .= rhs[isrt:iend, :]
    @. u += k1 * sixth_dt

    @. tmp = u_p + half_dt * k1
    level.t = t + half_dt
    f(level, rhs, tmp)
    k2[isrt:iend, :] .= rhs[isrt:iend, :]
    @. u += third_dt * k2

    @. tmp = u_p + half_dt * k2
    level.t = t + half_dt
    f(level, rhs, tmp)
    k3[isrt:iend, :] .= rhs[isrt:iend, :]
    @. u += third_dt * k3

    @. tmp = u_p + k3
    level.t = t + dt
    f(level, rhs, tmp)
    k4[isrt:iend, :] .= rhs[isrt:iend, :]
    @. u += sixth_dt * k4
    return level.t = t + dt
end

function rk4_dense_output_y!(y, theta, h, yn, k)
    theta2 = theta * theta
    theta3 = theta2 * theta
    ck1 = theta - 3 * theta2 / 2 + 2 * theta3 / 3
    ck23 = theta2 - 2 * theta3 / 3
    ck4 = -theta2 / 2 + 2 * theta3 / 3
    @. y = yn + (ck1 * k[1] + ck23 * (k[2] + k[3]) + ck4 * k[4]) * h
    return nothing
end

function rk4_dense_output_dy!(dy, theta, h, yn, k)
    theta2 = theta * theta
    ck1 = 1 - 3 * theta + 2 * theta2
    ck23 = 2 * theta - 2 * theta2
    ck4 = -theta + 2 * theta2
    @. dy = ck1 * k[1] + ck23 * (k[2] + k[3]) + ck4 * k[4]
    return nothing
end

function rk4_dense_output_d2y!(d2y, theta, h, yn, k)
    ck1 = -3 + 4 * theta
    ck23 = 2 - 4 * theta
    ck4 = -1 + 4 * theta
    coeff = 1 / h
    @. d2y = (ck1 * k[1] + ck23 * (k[2] + k[3]) + ck4 * k[4]) * coeff
    return nothing
end

function rk4_dense_output_d3y!(d3y, theta, h, yn, k)
    coeff = 4 / (h * h)
    @. d3y = (k[1] - k[2] - k[3] + k[4]) * coef
    return nothing
end
