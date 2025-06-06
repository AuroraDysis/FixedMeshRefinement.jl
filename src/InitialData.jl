module InitialData

include("Derivs.jl")
include("Sync.jl")
include("ODESolver.jl")
include("Physical.jl")

#===============================================================================
Initial Data Types:
    * Gaussian
===============================================================================#
function Gaussian!(gfs; amp = 1.0, sig = 0.25, x0 = 0.0)
    lmax = length(gfs.levels)
    for l = 1:lmax
        psi = gfs.levels[l].u[1]
        Pi = gfs.levels[l].u[2]
        x = gfs.levels[l].x
        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        Sync.restriction(gfs, l)
    end
end

function sinusoidal!(gfs)
    lmax = length(gfs.levels)
    for l = 1:lmax
        psi = gfs.levels[l].u[1]
        Pi = gfs.levels[l].u[2]
        x = gfs.levels[l].x
        @. psi = sin(2 * pi * (x - 0.0))
        @. Pi = -2 * pi * cos(2 * pi * (x - 0.0))
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        Sync.restriction(gfs, l)
    end
end

#===============================================================================
Spectial Treatment for prolongation
    * evolve backwards to file u_p
===============================================================================#
function NegativeWaveRHS!(level, r, u)
    Physical.WaveRHS!(level, r, u)
    @. r = -r
end

function MarchBackwards!(gfs)
    for l = 1:length(gfs.levels)
        if l > 1
            Sync.prolongation(gfs, l, false)
        end
        ODESolver.rk4!(NegativeWaveRHS!, gfs.levels[l])
        # save new u(-dt) -> u_p, u(0) -> u
        u = gfs.levels[l].u
        u_p = gfs.levels[l].u_p
        u_pp = gfs.levels[l].u_pp
        @. u_pp = u_p
        @. u_p = u
        @. u = u_pp
        gfs.levels[l].level.time = 0.0
    end
    gfs.grid.time = 0.0
end

end # module InitialData
