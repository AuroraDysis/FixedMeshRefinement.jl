function step!(
    grid::Grid{NumState,NumDiagnostic}; mongwane=false, apply_trans_zone=false
) where {NumState,NumDiagnostic}
    max_level = length(grid.levels)

    #-------------------------------------------------#
    # march the first substep for all levels          #
    #-------------------------------------------------#
    for l in 1:max_level  # notice that we march coarse level first
        if l > 1
            if mongwane
                prolongation_mongwane(grid, l, false)
            else
                prolongation(grid, l, false)
            end
            if apply_trans_zone
                apply_transition_zone(grid, l, false)
            end
        end
        mongwane ? rk4_mongwane!(f, grid.levels[l]) : rk4!(f, grid.levels[l])
    end

    #-------------------------------------------------#
    # march the other substeps to the same time slice #
    #-------------------------------------------------#
    if grid.subcycling
        levels = grid.levels
        dt_min = levels[max_level].dt
        substeps = ones(Int, max_level)
        for s in 2:(2^(max_level - 1))  # from second to final substep (of the finest level)
            for l in 2:max_level  # march all levels except the coarest (from coarse to fine)
                if l == max_level || (
                    isapprox(levels[l].time, levels[l + 1].time; rtol=1e-12) &&
                    abs(levels[l].time - levels[1].time) > dt_min
                )
                    substeps[l] += 1
                    if l < max_level
                        restriction(grid, l; apply_trans_zone=apply_trans_zone)  # from l+1 to l
                    end
                    # from l-1 to l
                    if mongwane
                        prolongation_mongwane(grid, l, mod(substeps[l], 2) == 0)
                    else
                        prolongation(grid, l, mod(substeps[l], 2) == 0)
                    end
                    if apply_trans_zone
                        apply_transition_zone(grid, l, mod(substeps[l], 2) == 0)
                    end
                    mongwane ? rk4_mongwane!(f, grid.levels[l]) : rk4!(f, grid.levels[l])
                end
            end
        end
    end

    #------------------------#
    # restriction all levels #
    #------------------------#
    for l in (max_level - 1):-1:1  # notice that we restrict fine level first
        restriction(grid, l; apply_trans_zone=apply_trans_zone)
    end

    #------------------#
    # update grid time #
    #------------------#
    grid.time = grid.levels[1].time

    return nothing
end

function rk4!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    u_pp = levf.u_pp
    r = levf.rhs
    w = levf.w
    lev = levf.lev
    t = lev.time
    dt = lev.dt

    @. u_pp = u_p * 1.0
    @. u_p = u * 1.0
    lev.time = t
    f(lev, r, u)
    @. u += r * (dt / 6)

    @. w = u_p + r * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    @. u += r * (dt / 3)

    @. w = u_p + r * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    @. u += r * (dt / 3)

    @. w = u_p + r * dt
    lev.time = t + dt
    f(lev, r, w)
    @. u += r * (dt / 6)
    lev.time = t + dt
end

function rk4_mongwane!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    r = levf.rhs
    w = levf.w
    k1 = levf.k[1]
    k2 = levf.k[2]
    k3 = levf.k[3]
    k4 = levf.k[4]
    lev = levf.lev
    t = lev.time
    dt = lev.dt
    isrt = lev.is_lev1 ? 1 : 1 + lev.nbuf
    iend = lev.is_lev1 ? lev.nxa : lev.nxa - lev.nbuf

    @. u_p = u * 1.0
    lev.time = t
    f(lev, r, u)
    for v = 1:levf.nd
        k1[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k1 * (1 / 6)

    @. w = u_p + k1 * (1 / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    for v = 1:levf.nd
        k2[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k2 * (1 / 3)

    @. w = u_p + k2 * (1 / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    for v = 1:levf.nd
        k3[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k3 * (1 / 3)

    @. w = u_p + k3
    lev.time = t + dt
    f(lev, r, w)
    for v = 1:levf.nd
        k4[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k4 * (1 / 6)
    lev.time = t + dt
end
