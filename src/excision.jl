export excise_level!, excise_grid!

"""
    excise_level!(level::Level, num_excise_points::NTuple{2,Float64})

Excise the boundary of a `Level` by a given number of points, as a prerequisite, the level must align with the physical boundary.
"""
function excise_level!(level::Level, num_excise_points::NTuple{2,Float64})
    (; is_physical_boundary) = level

    if (num_excise_points[1] > 0 && !is_physical_boundary[1]) ||
        (num_excise_points[2] > 0 && !is_physical_boundary[2])
        println("level.domain_box = ", level.domain_box)
        println("level.is_physical_boundary = ", level.is_physical_boundary)
        error("Level is not aligned with the physical boundary")
    end
end

"""
    move_physical_boundary!(grid::Grid, num_excise_points::NTuple{2,Float64})

Excise the grid by a given number of points; as a prerequisite, all levels must align with the physical boundary.
"""
function excise_grid!(grid::Grid, num_excise_points::NTuple{2,Float64})
    num_excise_points[1] >= 0 && num_excise_points[2] >= 0 ||
        error("num_excise_points must be non-negative")

    (; num_levels, levels) = grid

    # make sure all levels are at same time
    t = levels[1].t
    for l in 2:num_levels
        level = levels[l]
        isapprox_tol(level.t, t) ||
            error("Level $(l) is not at same time, t = $(t), levels[l].t = $(levels[l].t)")
    end

    for l in 1:num_levels
        level = levels[l]
        left_excise_points = num_excise_points[1] * 2^(l - 1)
        right_excise_points = num_excise_points[2] * 2^(l - 1)
        excise_level!(level, (left_excise_points, right_excise_points))
    end
end
