export rk4!

"""
    fill_buffer!(u, level::Level, stage::Int)

Fill the ghost cells of the solution array `u` from the `Yn_buffer` at a specific
Runge-Kutta `stage`. This is part of the subcycling-in-time method. After copying,
the buffer is filled with `NaN` to avoid accidental reuse.

# Arguments
- `u`: The solution array with ghost cells to be filled.
- `level::Level`: The grid level.
- `stage::Int`: The RK stage index.
- `fill_with_nan::Bool`: If `true`, fill the buffer with `NaN` to avoid accidental reuse.
  Defaults to `true`.
"""
function fill_buffer!(u, level::Level, stage::Int; fill_with_nan::Bool=true)
    (; Yn_buffer, num_boundary_points, is_physical_boundary) = level
    Yn = Yn_buffer[stage]
    boundary_indices = get_boundary_indices(level)
    for dir in 1:2
        if is_physical_boundary[dir]
            continue
        end
        for i in 1:num_boundary_points[dir]
            idx = boundary_indices[dir][i]
            u[:, idx] .= @view(Yn[:, i, dir])
        end
    end

    # fill with NaN to avoid accidental reuse
    if fill_with_nan
        Yn .= NaN
    end

    return nothing
end

"""
    rk4!(level::Level, f::Function, p; mongwane::Bool=false)

Perform a single step of the classic 4th-order Runge-Kutta (RK4) method on a given
`level`. It advances the solution by one time step `dt`.

# Arguments
- `level::Level`: The grid level to be updated.
- `f::Function`: The function that computes the right-hand side of the ODEs.
  It should have the signature `f(level, k, u, p, t)`.
- `p`: Parameters to be passed to the RHS function `f`.
- `mongwane::Bool`: If `true`, enables special buffer filling for Mongwane's
  subcycling method. Defaults to `false`.
"""
function rk4!(level::Level, f::Function, p; mongwane::Bool=false)
    (; t, dt) = level

    tmp = get_rk4_tmp_state(level)
    k1 = get_rk_stage(level, 1)
    k2 = get_rk_stage(level, 2)
    k3 = get_rk_stage(level, 3)
    k4 = get_rk_stage(level, 4)

    cycle_state!(level)

    u_p = get_state(level, -1)

    sixth_dt = dt / 6
    half_dt = dt / 2

    idx = get_rk4_evaluation_indices(level)

    f(level, k1, u_p, p, t)

    @.. tmp[:, idx] = u_p[:, idx] + half_dt * k1[:, idx]
    if mongwane
        fill_buffer!(tmp, level, 2)
    end
    f(level, k2, tmp, p, t + half_dt)

    @.. tmp[:, idx] = u_p[:, idx] + half_dt * k2[:, idx]
    if mongwane
        fill_buffer!(tmp, level, 3)
    end
    f(level, k3, tmp, p, t + half_dt)

    @.. tmp[:, idx] = u_p[:, idx] + dt * k3[:, idx]
    if mongwane
        fill_buffer!(tmp, level, 4)
    end
    f(level, k4, tmp, p, t + dt)

    u = get_state(level)
    @.. u[:, idx] = u_p[:, idx] + sixth_dt * (2 * (k2[:, idx] + k3[:, idx]) + (k1[:, idx] + k4[:, idx]))

    # update time
    level.t = t + dt

    return nothing
end
