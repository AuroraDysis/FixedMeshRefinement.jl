#===============================================================================
Initial Data Types:
    * Gaussian
===============================================================================#
function gaussian!(grid; amp=1.0, sig=0.25, x0=0.0)
    (; levels, num_levels) = grid
    for l in 1:num_levels
        psi = levels[l].u[1]
        Pi = levels[l].u[2]
        x = levels[l].x
        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end

    # restriction for consistence
    for l in (num_levels - 1):-1:1
        restriction!(grid, l)
    end
end

function sinusoidal!(grid)
    lmax = length(grid.levels)
    for l in 1:lmax
        psi = grid.levels[l].u[1]
        Pi = grid.levels[l].u[2]
        x = grid.levels[l].x
        @. psi = sin(2 * pi * (x - 0.0))
        @. Pi = -2 * pi * cos(2 * pi * (x - 0.0))
    end
    # restriction for consistence
    for l in (lmax - 1):-1:1
        restriction!(grid, l)
    end
end

#===============================================================================
Spectial Treatment for prolongation!
    * evolve backwards to file u_p
===============================================================================#
function NegativeWaveRHS!(level, r, u)
    Physical.rhs_wave!(level, r, u)
    @. r = -r
end

function MarchBackwards!(grid)
    for l in 1:length(grid.levels)
        if l > 1
            prolongation!(grid, l, false)
        end
        rk4!(grid.levels[l], NegativeWaveRHS!)
        # save new u(-dt) -> u_p, u(0) -> u
        u = grid.levels[l].u
        u_p = grid.levels[l].u_p
        u_pp = grid.levels[l].u_pp
        @. u_pp = u_p
        @. u_p = u
        @. u = u_pp
        grid.levels[l].t = 0.0
    end
    return grid.t = 0.0
end
