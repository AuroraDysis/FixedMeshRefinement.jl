export integrate_trapezoidal

"""
    integrate_trapezoidal(grid::Grid, getter::Function)

Integrate a variable over all levels of the grid using the trapezoidal rule.
The `getter` function specifies which data to retrieve the data of the variable from each level.
`getter(level)` should return an `OffsetArray` from which data is read.
"""
function integrate_trapezoidal(grid::Grid, getter::Function)
    x, y = merge_grid_levels(grid, getter)

    n = length(x)
    length(x) == length(y) || error("x and y vectors must be of the same length!")
    n ≥ 2 || error("vectors must contain at least two elements")

    retval = 0.0
    for i in 1:blocklength(x)
        block_x = getblock(x, i)
        block_y = getblock(y, i)

        dx = block_x[begin+1] - block_x[begin]

        y_prev = i == 1 ? block_y[1] : getblock(y, i - 1)[end]
        local_retval = i == 1 ? block_y[2] : block_y[1]
        start_idx = i == 1 ? 3 : 2
        @inbounds @fastmath @simd for j in start_idx:(length(block_y) - 1)
            local_retval += block_y[j]
        end
        @inbounds retval += dx * (local_retval + (y_prev + block_y[end]) / 2)
    end

    return retval
end
