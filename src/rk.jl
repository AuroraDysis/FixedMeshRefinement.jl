function rk4!(level::Level, f::Function)
    (; u, u_p, tmp, k, t, dt) = level
    sixth_dt = dt / 6
    half_dt = dt / 2

    k1, k2, k3, k4 = k

    f(level, k1, u, t)

    @. tmp = u_p + half_dt * k1
    f(level, k2, tmp, t + half_dt)

    @. tmp = u_p + half_dt * k2
    f(level, k3, tmp, t + half_dt)

    @. tmp = u_p + dt * k3
    f(level, k4, tmp, t + dt)

    @. u = u_p + sixth_dt * (2 * (k2 + k3) + (k1 + k4))

    # update time
    level.t = t + dt

    return nothing
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
