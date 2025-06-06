#===============================================================================
Functions needed by Mongwane's subcycling method
===============================================================================#
function calc_kfs_from_kcs(kcs, dtc, interp_in_time::Bool)
    t0_f = interp_in_time ? 0.5 : 0.0
    dtf = 0.5 * dtc
    d1yc = DenseOutput.dy1(t0_f, dtc, kcs)
    d2yc = DenseOutput.dy2(t0_f, dtc, kcs)
    d3yc = DenseOutput.dy3(t0_f, dtc, kcs)
    fyd2yc = 4 * (kcs[3] - kcs[2]) / dtc^3
    return [
        dtf * d1yc,
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc - fyd2yc),
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc + fyd2yc),
    ]
end

function transition_profile(xl, xh, x; type = 1)
    t0 = (x - xl) / (xh - xl)
    t = t0 < 0 ? 0 : (t0 > 1 ? 1 : t0)

    if type == 1  # boxstep
        return t
    elseif type == 2  # smoothstep
        return 3 * t^2 - 2 * t^3
    elseif type == 3  # smootherstep
        return 10 * t^3 - 15 * t^4 + 6 * t^5
    else
        println("Transition profile (type $type) is not supported yet.")
        exit()
    end
end

#===============================================================================
apply_transition_zone: apply transition zone
===============================================================================#
function apply_transition_zone(grid, l, interp_in_time::Bool)
    num_total_points = grid.levels[l].num_total_points
    num_buffer_points = grid.levels[l].num_buffer_points
    num_transition_points = grid.levels[l].num_transition_points
    spatial_interpolation_order = grid.levels[l].spatial_interpolation_order
    parent_map = grid.levels[l].parent_map
    is_aligned = grid.levels[l].is_aligned
    levf = grid.levels[l]
    levc = grid.levels[l-1]
    # for transition zone
    domain_box = grid.levels[l].domain_box
    dxf = grid.levels[l].dx
    @assert(isapprox(domain_box[1], grid.levels[l].x[1+num_buffer_points]; rtol = 1e-12))
    @assert(isapprox(domain_box[2], grid.levels[l].x[num_total_points-num_buffer_points]; rtol = 1e-12))

    for j = 1:2  # left or right
        a = (j == 1) ? domain_box[1] : domain_box[2]
        b = (j == 1) ? domain_box[1] + (num_transition_points - 1) * dxf : domain_box[2] - (num_transition_points - 1) * dxf
        for v = 1:grid.NumState
            uf = levf.state[v]
            uc_p = levc.state_prev[v]
            for i = 1:num_transition_points
                f = (j == 1) ? i + num_buffer_points : num_total_points - i + 1 - num_buffer_points
                c = parent_map[f]
                w = transition_profile(a, b, grid.levels[l].x[f])
                if is_aligned[f]
                    kcs = [levc.k[m][v][c] for m = 1:4]
                    ys = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                    uf[f] = (1 - w) * ys + w * uf[f]
                else
                    nys = spatial_interpolation_order + 1
                    ioffset = (mod(nys, 2) == 0) ? div(nys, 2) : div(nys, 2) + 1
                    ys = zeros(Float64, nys)
                    for ic = 1:nys
                        ic_grid = c + ic - ioffset;
                        kcs = [levc.k[m][v][ic_grid] for m = 1:4]
                        ys[ic] =
                            interp_in_time ? DenseOutput.y(0.5, uc_p[ic_grid], kcs) :
                            uc_p[ic_grid]
                    end
                    uf[f] = (1 - w) * Algo.Interpolation(ys, ioffset, spatial_interpolation_order) + w * uf[f]
                end
            end
        end
    end
end

#===============================================================================
prolongation_mongwane: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function prolongation_mongwane(grid, l, interp_in_time::Bool)
    num_total_points = grid.levels[l].num_total_points
    num_buffer_points = grid.levels[l].num_buffer_points
    spatial_interpolation_order = grid.levels[l].spatial_interpolation_order
    parent_map = grid.levels[l].parent_map
    is_aligned = grid.levels[l].is_aligned
    dtc = grid.levels[l-1].dt
    levf = grid.levels[l]
    levc = grid.levels[l-1]

    for j = 1:2  # left or right
        for v = 1:grid.NumState
            uf = levf.state[v]
            uc_p = levc.state_prev[v]
            for i = 1:num_buffer_points
                f = (j == 1) ? i : num_total_points - i + 1
                c = parent_map[f]
                if is_aligned[f]
                    kcs = [levc.k[m][v][c] for m = 1:4]
                    kfs = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                    # setting k
                    for m = 1:3
                        levf.k[m][v][f] = kfs[m]
                    end
                    # setting state
                    uf[f] = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                else
                    nys = spatial_interpolation_order + 1
                    ioffset = (mod(nys, 2) == 0) ? div(nys, 2) : div(nys, 2) + 1
                    kfss = zeros(Float64, 3, nys)
                    ys = zeros(Float64, nys)
                    for ic = 1:nys
                        ic_grid = c + ic - ioffset;
                        kcs = [levc.k[m][v][ic_grid] for m = 1:4]
                        kfss[:, ic] = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                        ys[ic] =
                            interp_in_time ? DenseOutput.y(0.5, uc_p[ic_grid], kcs) :
                            uc_p[ic_grid]
                    end
                    # setting k
                    for m = 1:3
                        levf.k[m][v][f] = Algo.Interpolation(kfss[m, :], ioffset, spatial_interpolation_order)
                    end
                    # setting state
                    uf[f] = Algo.Interpolation(ys, ioffset, spatial_interpolation_order)
                end
            end
        end
    end
end

#===============================================================================
prolongation:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function prolongation(grid, l, interp_in_time::Bool; ord_t = 2)
    num_total_points = grid.levels[l].num_total_points
    num_buffer_points = grid.levels[l].num_buffer_points
    spatial_interpolation_order = grid.levels[l].spatial_interpolation_order
    parent_map = grid.levels[l].parent_map
    is_aligned = grid.levels[l].is_aligned
    levf = grid.levels[l]
    levc = grid.levels[l-1]

    for j = 1:2  # left or right
        for v = 1:grid.NumState
            uf = levf.state[v]
            uc_p = levc.state_prev[v]
            if interp_in_time
                uc = levc.state[v]
                uc_pp = levc.state_prev_prev[v]
                for i = 1:num_buffer_points
                    f = (j == 1) ? i : num_total_points - i + 1
                    c = parent_map[f]
                    if is_aligned[f]
                        uf[f] = Algo.Interpolation([uc_pp[c], uc_p[c], uc[c]], 2, ord_t)
                    else
                        nucss = spatial_interpolation_order + 1
                        ioffset = (mod(nucss, 2) == 0) ? div(nucss, 2) : div(nucss, 2) + 1
                        ucss = zeros(Float64, 3, nucss)
                        for ic = 1:nucss
                            ic_grid = c + ic - ioffset;
                            ucss[:, ic] = [uc_pp[ic_grid], uc_p[ic_grid], uc[ic_grid]]
                        end
                        uf[f] = Algo.Interpolation(
                            [Algo.Interpolation(ucss[m, :], ioffset, spatial_interpolation_order) for m = 1:3],
                            2,
                            ord_t,
                        )
                    end
                end
            else
                for i = 1:num_buffer_points
                    f = (j == 1) ? i : num_total_points - i + 1
                    c = parent_map[f]
                    uf[f] = ((is_aligned[f]) ? uc_p[c] : Algo.Interpolation(uc_p, c, spatial_interpolation_order))
                end
            end
        end
    end
end

#===============================================================================
restriction:
    * from level l+1 to level l
    * we assume that we always march fine level first (for l in lmax-1:-1:1)
    * we assume all the levels are at the same time slice
===============================================================================#
function restriction(grid, l; apply_trans_zone = false)
    num_total_points = grid.levels[l+1].num_total_points
    num_buffer_points = grid.levels[l+1].num_buffer_points
    num_transition_points = grid.levels[l+1].num_transition_points
    parent_map = grid.levels[l+1].parent_map
    is_aligned = grid.levels[l+1].is_aligned
    isrt = apply_trans_zone ? 1 + num_buffer_points + num_transition_points : 1 + num_buffer_points
    iend = apply_trans_zone ? num_total_points - num_buffer_points - num_transition_points : num_total_points - num_buffer_points
    for v = 1:grid.NumState
        uf = grid.levels[l+1].state[v]
        uc = grid.levels[l].state[v]
        for f = isrt:iend  # only interior
            if is_aligned[f]
                uc[parent_map[f]] = uf[f]
            end
        end
    end
end