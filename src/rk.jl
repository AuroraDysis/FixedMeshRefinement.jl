function fill_buffer!(u, level::Level, stage::Int)
    (; Yn_buffer, buffer_indices, num_buffer_points, x, physical_domain_box) = level
    Yn = Yn_buffer[stage]
    for dir in 1:2, i in 1:num_buffer_points
        idx = buffer_indices[dir][i]
        x_pos = x[idx]
        # don't change if the points are outside the physical boundary
        if x_pos < physical_domain_box[1] || x_pos > physical_domain_box[2]
            continue
        end
        u[idx, :] .= @view(Yn[idx, :, dir])
    end
    fill!(Yn, NaN)
    return nothing
end

function rk4!(level::Level, f::Function; mongwane::Bool=false)
    (; state, tmp, k, t, dt) = level

    cycle_state!(level)

    u_p = state[end - 1]

    sixth_dt = dt / 6
    half_dt = dt / 2

    k1, k2, k3, k4 = k

    f(level, k1, u_p, t)

    @. tmp = u_p + half_dt * k1
    if mongwane
        fill_buffer!(tmp, level, 2)
    end
    f(level, k2, tmp, t + half_dt)

    @. tmp = u_p + half_dt * k2
    if mongwane
        fill_buffer!(tmp, level, 3)
    end
    f(level, k3, tmp, t + half_dt)

    @. tmp = u_p + dt * k3
    if mongwane
        fill_buffer!(tmp, level, 4)
    end
    f(level, k4, tmp, t + dt)

    u = state[end]
    @. u = u_p + sixth_dt * (2 * (k2 + k3) + (k1 + k4))

    # update time
    level.t = t + dt

    return nothing
end
