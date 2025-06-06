module Algo

function Interpolation(state, i, ord)
    if ord == 1
        return (state[i] + state[i+1]) * 0.5
    elseif ord == 2
        return (-state[i-1] + 6 * state[i] + 3 * state[i+1]) * 0.125
    elseif ord == 3
        return (-state[i-1] + 9 * state[i] + 9 * state[i+1] - state[i+2]) * 0.0625
    elseif ord == 5
        return (
            3 * state[i-2] - 25 * state[i-1] + 150 * state[i] + 150 * state[i+1] - 25 * state[i+2] + 3 * state[i+3]
        ) / 256
    else
        println("Interpolation order not supported yet: ord = ", ord)
        exit()
    end
end

end
