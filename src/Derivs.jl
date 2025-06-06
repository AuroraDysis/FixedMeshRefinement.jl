module Derivs

#===============================================================================
Derivatives:
    * single point
===============================================================================#
function deriv_1st(state, i, dx, ord)
    if ord == 2
        return (-state[i-1] + state[i+1]) / (2 * dx)
    elseif ord == 4
        return (state[i-2] - 8 * state[i-1] + 8 * state[i+1] - state[i+2]) / (12 * dx)
    else
        println("Finite difference order not supported yet: ord = ", ord)
        exit()
    end
end

function deriv_2nd(state, i, dx, ord)
    if ord == 2
        return (state[i-1] - 2 * state[i] + state[i+1]) / (dx * dx)
    elseif ord == 4
        return (-state[i-2] + 16 * state[i-1] - 30 * state[i] + 16 * state[i+1] - state[i+2]) / (12 * dx * dx)
    else
        println("Finite difference order not supported yet: ord = ", ord)
        exit()
    end
end

function deriv_diss(state, i, dx, diss_ord)
    if diss_ord == 4
        return ((state[i+2] + state[i-2]) - 4 * (state[i+1] + state[i-1]) + 6 * state[i]) / dx
    elseif diss_ord == 6
        return (
            (state[i+3] + state[i-3]) - 6 * (state[i+2] + state[i-2]) + 15 * (state[i+1] + state[i-1]) - 20 * state[i]
        ) / dx
    else
        println("KO dissipation order not supported yet: ord = ", ord)
        exit()
    end
end

#===============================================================================
Derivatives:
    * single level
===============================================================================#
function derivs_1st!(du, state, dx, ord)
    istart = 1 + div(ord, 2)
    iend = length(state) - div(ord, 2)
    for i = istart:iend
        du[i] = deriv_1st(state, i, dx, ord)
    end
end

function derivs_2nd!(ddu, state, dx, ord)
    istart = 1 + div(ord, 2)
    iend = length(state) - div(ord, 2)
    for i = istart:iend
        ddu[i] = deriv_2nd(state, i, dx, ord)
    end
end

function derivs_diss!(dissipation, state, dx, ord)
    diss_ord = ord + 2
    sign = (mod(diss_ord, 4) == 0 ? -1 : +1)
    istart = 1 + div(diss_ord, 2)
    iend = length(state) - div(diss_ord, 2)
    for i = istart:iend
        dissipation[i] = deriv_diss(state, i, dx, diss_ord) * (sign / 2^(diss_ord))
    end
end

#===============================================================================
Derivatives:
    * single level
    * for varlist
===============================================================================#
function calc_du!(du::Array{Array{Float64,1},1}, state::Array{Array{Float64,1},1}, dx, ord)
    for i = 1:length(state)
        derivs_1st!(du[i], state[i], dx, ord)
    end
end

function calc_ddu!(ddu::Array{Array{Float64,1},1}, state::Array{Array{Float64,1},1}, dx, ord)
    for i = 1:length(state)
        derivs_2nd!(ddu[i], state[i], dx, ord)
    end
end

end # module Derivs
