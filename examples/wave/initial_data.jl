#===============================================================================
Initial Data Types:
    * Gaussian
===============================================================================#
function gaussian!(grid; amp=1.0, sig=0.25, x0=0.0)
    (; levels, num_levels) = grid
    for l in 1:num_levels
        level = levels[l]
        u = get_state(level)
        psi = @view(u[:, 1])
        Pi = @view(u[:, 2])
        x = get_x(level)
        @.. psi = amp * exp(-((x - x0) / sig)^2)
        @.. Pi = 0.0
    end

    # restriction for consistence
    for l in (num_levels - 1):-1:1
        restrict_injection!(grid, l)
    end
end

function sinusoidal!(grid)
    (; levels, num_levels) = grid

    for l in 1:num_levels
        level = levels[l]

        u = get_state(level)
        psi = @view(u[:, 1])
        Pi = @view(u[:, 2])

        x = get_x(level)
        @.. psi = sin(2 * pi * (x - 0.0))
        @.. Pi = -2 * pi * cos(2 * pi * (x - 0.0))
    end

    # restriction for consistence
    for l in (num_levels - 1):-1:1
        restrict_injection!(grid, l)
    end
end

#===============================================================================
Spectial Treatment for prolongate!
    * evolve backwards to file u_p
===============================================================================#
function wave_rhs_backward!(level, r, u, p, t)
    wave_rhs!(level, r, u, p, t)
    @.. r = -r
end

function march_backwards!(grid, p)
    (; levels, num_levels) = grid

    for l in 1:num_levels
        if l > 1
            prolongate!(grid, l, false)
        end
        level = levels[l]
        rk4!(level, wave_rhs_backward!, p)

        # save new u(-dt) -> u_p, u(0) -> u
        tmp = get_rk4_tmp_state(level)
        u = get_state(level)
        u_p = get_state(level, -1)
        tmp .= u
        u .= u_p
        u_p .= tmp

        level.t = 0.0
    end

    # reset time
    grid.t = 0.0
    for l in 1:num_levels
        levels[l].t = 0.0
    end

    return nothing
end
