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

function transition_profile(xl, xh, x; type=1)
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
apply_transition_zone!: apply transition zone
===============================================================================#
function apply_transition_zone!(grid, l, interp_in_time::Bool)
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_total_points,
        num_buffer_points,
        num_transition_points,
        spatial_interpolation_order,
        parent_indices,
    ) = fine_level

    # for transition zone
    domain_box = fine_level.domain_box
    dxf = grid.levels[l].dx
    @assert(isapprox(domain_box[1], grid.levels[l].x[1 + num_buffer_points]; rtol=1e-12))
    @assert(
        isapprox(
            domain_box[2],
            grid.levels[l].x[num_total_points - num_buffer_points];
            rtol=1e-12,
        )
    )

    for j in 1:2  # left or right
        a = (j == 1) ? domain_box[1] : domain_box[2]
        b = if (j == 1)
            domain_box[1] + (num_transition_points - 1) * dxf
        else
            domain_box[2] - (num_transition_points - 1) * dxf
        end
        for v in 1:(grid.NumState)
            uf = fine_level.u[v]
            uc_p = coarse_level.u_p[v]
            for i in 1:num_transition_points
                f = if (j == 1)
                    i + num_buffer_points
                else
                    num_total_points - i + 1 - num_buffer_points
                end
                c = parent_map[f]
                w = transition_profile(a, b, grid.levels[l].x[f])
                if is_aligned[f]
                    kcs = [coarse_level.k[m][v][c] for m in 1:4]
                    ys = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                    uf[f] = (1 - w) * ys + w * uf[f]
                else
                    nys = spatial_interpolation_order + 1
                    ioffset = (mod(nys, 2) == 0) ? div(nys, 2) : div(nys, 2) + 1
                    ys = zeros(Float64, nys)
                    for ic in 1:nys
                        ic_grid = c + ic - ioffset
                        kcs = [coarse_level.k[m][v][ic_grid] for m in 1:4]
                        ys[ic] = if interp_in_time
                            DenseOutput.y(0.5, uc_p[ic_grid], kcs)
                        else
                            uc_p[ic_grid]
                        end
                    end
                    uf[f] =
                        (1 - w) *
                        Algo.Interpolation(ys, ioffset, spatial_interpolation_order) +
                        w * uf[f]
                end
            end
        end
    end
end

#===============================================================================
prolongation_mongwane!: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function prolongation_mongwane!(grid, l, interp_in_time::Bool)
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (; num_total_points, num_buffer_points, spatial_interpolation_order, parent_indices) =
        fine_level

    dtc = coarse_level.dt

    for j in 1:2  # left or right
        for v in 1:(grid.NumState)
            uf = fine_level.u[v]
            uc_p = coarse_level.u_p[v]
            for i in 1:num_buffer_points
                f = (j == 1) ? i : num_total_points - i + 1
                c = parent_map[f]
                if is_aligned[f]
                    kcs = [coarse_level.k[m][v][c] for m in 1:4]
                    kfs = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                    # setting k
                    for m in 1:3
                        fine_level.k[m][v][f] = kfs[m]
                    end
                    # setting u
                    uf[f] = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                else
                    nys = spatial_interpolation_order + 1
                    ioffset = (mod(nys, 2) == 0) ? div(nys, 2) : div(nys, 2) + 1
                    kfss = zeros(Float64, 3, nys)
                    ys = zeros(Float64, nys)
                    for ic in 1:nys
                        ic_grid = c + ic - ioffset
                        kcs = [coarse_level.k[m][v][ic_grid] for m in 1:4]
                        kfss[:, ic] = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                        ys[ic] = if interp_in_time
                            DenseOutput.y(0.5, uc_p[ic_grid], kcs)
                        else
                            uc_p[ic_grid]
                        end
                    end
                    # setting k
                    for m in 1:3
                        fine_level.k[m][v][f] = Algo.Interpolation(
                            kfss[m, :], ioffset, spatial_interpolation_order
                        )
                    end
                    # setting u
                    uf[f] = Algo.Interpolation(ys, ioffset, spatial_interpolation_order)
                end
            end
        end
    end
end

#===============================================================================
prolongation!:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function prolongation!(grid::Grid, l::Int, interp_in_time::Bool; ord_t=2)
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (; num_total_points, num_buffer_points, spatial_interpolation_order, parent_indices) =
        fine_level

    # j: 1: left, 2: right
    for j in 1:2
        uf = fine_level.u[v]
        uc_p = coarse_level.u_p[v]
        if interp_in_time
            uc = coarse_level.u[v]
            uc_pp = coarse_level.u_pp[v]
            for i in 1:num_buffer_points
                f = (j == 1) ? i : num_total_points - i + 1
                c = parent_map[f]
                if is_aligned[f]
                    uf[f] = Algo.Interpolation([uc_pp[c], uc_p[c], uc[c]], 2, ord_t)
                else
                    nucss = spatial_interpolation_order + 1
                    ioffset = (mod(nucss, 2) == 0) ? div(nucss, 2) : div(nucss, 2) + 1
                    ucss = zeros(Float64, 3, nucss)
                    for ic in 1:nucss
                        ic_grid = c + ic - ioffset
                        ucss[:, ic] = [uc_pp[ic_grid], uc_p[ic_grid], uc[ic_grid]]
                    end
                    uf[f] = Algo.Interpolation(
                        [
                            Algo.Interpolation(
                                ucss[m, :], ioffset, spatial_interpolation_order
                            ) for m in 1:3
                        ],
                        2,
                        ord_t,
                    )
                end
            end
        else
            for i in 1:num_buffer_points
                f = (j == 1) ? i : num_total_points - i + 1
                c = parent_map[f]
                uf[f] = (
                    if (is_aligned[f])
                        uc_p[c]
                    else
                        Algo.Interpolation(uc_p, c, spatial_interpolation_order)
                    end
                )
            end
        end
    end
end
