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
        psi = grid.levels[l].state[1]
        Pi = grid.levels[l].state[2]
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
        psi = grid.levels[l].state[1]
        Pi = grid.levels[l].state[2]
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
Spectial Treatment for prolongation
    * evolve backwards to file state_prev
===============================================================================#
function NegativeWaveRHS!(level, r, state)
    Physical.WaveRHS!(level, r, state)
    @. r = -r
end

function MarchBackwards!(grid)
    for l = 1:length(grid.levels)
        if l > 1
            prolongation(grid, l, false)
        end
        ODESolver.rk4!(NegativeWaveRHS!, grid.levels[l])
        # save new state(-dt) -> state_prev, state(0) -> state
        state = grid.levels[l].state
        state_prev = grid.levels[l].state_prev
        state_prev_prev = grid.levels[l].state_prev_prev
        @. state_prev_prev = state_prev
        @. state_prev = state
        @. state = state_prev_prev
        grid.levels[l].level.time = 0.0
    end
    grid.time = 0.0
end

end # module InitialData
