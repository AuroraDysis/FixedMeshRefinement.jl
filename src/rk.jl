function rk4!(level::Level, f::Function; mongwane::Bool=false)
    (; state, tmp, k, t, dt, Yn_buffer, num_buffer_points, buffer_indices) = level

    cycle_state!(level)

    u = state[end]
    u_p = state[end - 1]

    sixth_dt = dt / 6
    half_dt = dt / 2

    k1, k2, k3, k4 = k

    if mongwane
        Y1 = Yn_buffer[1]
        for dir in 1:2
            u_p[buffer_indices[dir], :] .= @view(Y1[:, :, dir])
        end
    end
    f(level, k1, u_p, t)

    @. tmp = u_p + half_dt * k1
    if mongwane
        Y2 = Yn_buffer[2]
        for dir in 1:2
            tmp[buffer_indices[dir], :] .= @view(Y2[:, :, dir])
        end
    end
    f(level, k2, tmp, t + half_dt)

    @. tmp = u_p + half_dt * k2
    if mongwane
        Y3 = Yn_buffer[3]
        for dir in 1:2
            tmp[buffer_indices[dir], :] .= @view(Y3[:, :, dir])
        end
    end
    f(level, k3, tmp, t + half_dt)

    @. tmp = u_p + dt * k3
    if mongwane
        Y4 = Yn_buffer[4]
        for dir in 1:2
            tmp[buffer_indices[dir], :] .= @view(Y4[:, :, dir])
        end
    end
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
    @. d3y = (k[1] - k[2] - k[3] + k[4]) * coeff
    return nothing
end
