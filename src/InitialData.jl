module InitialData

include("Derivs.jl")
include("ODESolver.jl")
include("Physical.jl")

#===============================================================================
Initial Data Types:
    * Gaussian
===============================================================================#
function Gaussian!(grid; amp = 1.0, sig = 0.25, x0 = 0.0)
    lmax = length(grid.levels)
    for l = 1:lmax
        psi = grid.levels[l].u[1]
        Pi = grid.levels[l].u[2]
        x = grid.levels[l].x
        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        restriction!(grid, l)
    end
end

function sinusoidal!(grid)
    lmax = length(grid.levels)
    for l = 1:lmax
        psi = grid.levels[l].u[1]
        Pi = grid.levels[l].u[2]
        x = grid.levels[l].x
        @. psi = sin(2 * pi * (x - 0.0))
        @. Pi = -2 * pi * cos(2 * pi * (x - 0.0))
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        restriction!(grid, l)
    end
end

#===============================================================================
Spectial Treatment for prolongation!
    * evolve backwards to file u_p
===============================================================================#
function NegativeWaveRHS!(level, r, u)
    Physical.WaveRHS!(level, r, u)
    @. r = -r
end

function MarchBackwards!(grid)
    for l = 1:length(grid.levels)
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
    grid.t = 0.0
end

end # module InitialData
