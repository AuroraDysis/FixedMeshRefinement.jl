export integrate_trapezoidal, integrate_simpson, integrate_block_trapezoidal, integrate_block_simpson

function integrate_block_simpson(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}
    length(f) >= 4 || throw(ArgumentError("f must have at least 4 elements"))

    @inbounds val =
        (
            17 * (f[1] + f[end]) +
            59 * (f[2] + f[end - 1]) +
            43 * (f[3] + f[end - 2]) +
            49 * (f[4] + f[end - 3])
        ) / 48
    @inbounds @simd for i in 5:(length(f) - 4)
        val += f[i]
    end
    @inbounds return val * dx
end

function integrate_block_trapezoidal(f::AbstractVector{T}, dx::T) where {T<:AbstractFloat}
    length(f) >= 3 || throw(ArgumentError("f must have at least 3 elements"))

    @inbounds val = f[2]
    @inbounds @simd for i in 3:(length(f) - 1)
        val += f[i]
    end
    @inbounds return dx * (val + (f[1] + f[end]) / 2)
end

"""
    integrate_trapezoidal(grid::Grid, getter::Function)

Integrate a variable over all levels of the grid using the trapezoidal rule.
The `getter` function specifies which data to retrieve the data of the variable from each level.
`getter(level)` should return an `SubArray` of `OffsetArray`.
"""
function integrate_trapezoidal(grid::Grid, getter::Function)
    var_x, var_y = merge_grid_levels(grid, getter; include_overlap_points=true)

    retval = 0.0
    for i in 1:blocklength(var_x)
        x = getblock(var_x, i)
        y = getblock(var_y, i)

        dx = step(x)
        retval += integrate_block_trapezoidal(y, dx)
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
    var_x, var_y = merge_grid_levels(grid, getter; include_overlap_points=true)

    retval = 0.0
    for i in 1:blocklength(var_x)
        x = getblock(var_x, i)
        y = getblock(var_y, i)

        dx = step(x)
        retval += integrate_block_simpson(y, dx)
    end

    return retval
end
