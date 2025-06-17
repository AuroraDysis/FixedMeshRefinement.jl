export integrate_trapezoidal, integrate_simpson

"""
    integrate_trapezoidal(grid::Grid, getter::Function)

Integrate a variable over all levels of the grid using the trapezoidal rule.
The `getter` function specifies which data to retrieve the data of the variable from each level.
`getter(level)` should return an `SubArray` of `OffsetArray`.
"""
function integrate_trapezoidal(grid::Grid, getter::Function)
    var_x, var_y = merge_grid_levels(grid, getter)

    retval = 0.0
    for i in 1:blocklength(var_x)
        x = getblock(var_x, i)
        y = getblock(var_y, i)

        dx = x[begin + 1] - x[begin]

        first_idx = firstindex(y)
        last_idx = lastindex(y)

        y_prev = i == 1 ? y[first_idx] : getblock(var_y, i - 1)[end]
        local_retval = i == 1 ? y[first_idx + 1] : y[first_idx]
        start_idx = i == 1 ? first_idx + 2 : first_idx + 1
        @inbounds @fastmath @simd for j in start_idx:(last_idx - 1)
            local_retval += y[j]
        end
        @inbounds retval += dx * (local_retval + (y_prev + y[end]) / 2)
    end

    return retval
end

"""
    integrate_simpson(grid::Grid, getter::Function)

Integrate a variable over all levels of the grid using Simpson's rule.
The `getter` function specifies which data to retrieve the data of the variable from each level.
`getter(level)` should return an `SubArray` of `OffsetArray`.
"""
function integrate_simpson(grid::Grid, getter::Function)
    var_x, var_y = merge_grid_levels(grid, getter)

    retval = 0.0
    for i in 1:blocklength(var_x)
        x = getblock(var_x, i)
        y = getblock(var_y, i)

        dx = x[begin + 1] - x[begin]

        y_prev = i == 1 ? y[begin] : getblock(var_y, i - 1)[end]
        p1 = i == 1 ? y[begin+1] : y[begin]
        p2 = i == 1 ? y[begin+2] : y[begin+1]
        p3 = i == 1 ? y[begin+3] : y[begin+2]

        @inbounds local_retval = (
            17 * (y_prev + y[end]) +
            59 * (p1 + y[end-1]) +
            43 * (p2 + y[end-2]) +
            49 * (p3 + y[end-3])
        ) / 48

        first_idx = firstindex(y)
        last_idx = lastindex(y)
        start_idx = i == 1 ? first_idx + 4 : first_idx + 3
        @inbounds @fastmath @simd for j in start_idx:(last_idx-4)
            local_retval += y[j]
        end
        retval += dx * local_retval
    end

    return retval
end
