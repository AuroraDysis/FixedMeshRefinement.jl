function rk4!(f::Function, level)
    cycle_state!(level)

    u = level.u
    u_p = level.u_p
    u_pp = level.u_pp

    rhs = level.rhs
    tmp = level.tmp
    t = level.t
    dt = level.dt

    @. u_pp = u_p
    @. u_p = u
    level.t = t
    f(level, rhs, u)
    @. u += rhs * (dt / 6)

    @. tmp = u_p + rhs * (dt / 2)
    level.t = t + 0.5 * dt
    f(level, rhs, tmp)
    @. u += rhs * (dt / 3)

    @. tmp = u_p + rhs * (dt / 2)
    level.t = t + 0.5 * dt
    f(level, rhs, tmp)
    @. u += rhs * (dt / 3)

    @. tmp = u_p + rhs * dt
    level.t = t + dt
    f(level, rhs, tmp)
    @. u += rhs * (dt / 6)
    return level.t = t + dt
end

function rk4_mongwane!(
    f::Function, level::Level{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    (; u, u_p, rhs, tmp, k, t, dt) = level
    u = level.u
    u_p = level.u_p
    rhs = level.rhs
    tmp = level.tmp
    k1 = level.k[1]
    k2 = level.k[2]
    k3 = level.k[3]
    k4 = level.k[4]
    t = level.t
    dt = level.dt
    isrt = level.is_base_level ? 1 : 1 + level.num_buffer_points
    iend = if level.is_base_level
        level.num_total_points
    else
        level.num_total_points - level.num_buffer_points
    end

    @. u_p = u
    level.t = t
    f(level, rhs, u)
    for v in 1:NumState
        k1[v][isrt:iend] = rhs[v][isrt:iend] * dt
    end
    @. u += k1 * (1 / 6)

    @. tmp = u_p + k1 * (1 / 2)
    level.t = t + 0.5 * dt
    f(level, rhs, tmp)
    for v in 1:NumState
        k2[v][isrt:iend] = rhs[v][isrt:iend] * dt
    end
    @. u += k2 * (1 / 3)

    @. tmp = u_p + k2 * (1 / 2)
    level.t = t + 0.5 * dt
    f(level, rhs, tmp)
    for v in 1:NumState
        k3[v][isrt:iend] = rhs[v][isrt:iend] * dt
    end
    @. u += k3 * (1 / 3)

    @. tmp = u_p + k3
    level.t = t + dt
    f(level, rhs, tmp)
    for v in 1:NumState
        k4[v][isrt:iend] = rhs[v][isrt:iend] * dt
    end
    @. u += k4 * (1 / 6)
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
