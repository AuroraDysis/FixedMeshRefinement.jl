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

    @inbounds retval = (x[2] - x[1]) * (y[1] + y[2])
    @inbounds @fastmath @simd for i in 2:(n - 1)
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end

    return retval / 2
end
