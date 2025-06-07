function prolongation_spatial_interpolate!(res, u, i, order)
    if order == 1
        # {0.5, 0.5}
        @. res = (u[i] + u[i + 1]) * 0.5
    elseif order == 2
        # {-0.125, 0.75, 0.375}
        @. res = -0.125 * u[i - 1] + 0.75 * u[i] + 0.375 * u[i + 1]
    elseif order == 3
        # {-0.0625, 0.5625, 0.5625, -0.0625}
        @. res = -0.0625 * (u[i - 1] + u[i + 2]) + 0.5625 * (u[i] + u[i + 1])
    elseif order == 4
        # {0.0234375, -0.15625, 0.703125, 0.46875, -0.0390625}
        @. res =
            0.0234375 * u[i - 2] - 0.15625 * u[i - 1] +
            0.703125 * u[i] +
            0.46875 * u[i + 1] - 0.0390625 * u[i + 2]
    elseif order == 5
        # {0.0117188, -0.0976563, 0.585938, 0.585938, -0.0976563, 0.0117188}
        @. res =
            0.0117188 * (u[i - 2] + u[i + 3]) - 0.0976563 * (u[i - 1] + u[i + 2]) +
            0.585938 * (u[i] + u[i + 1])
    else
        println("Spatial interpolation order not supported yet: order = ", order)
        exit()
    end
end

# u from oldest to newest
function prolongation_time_interpolate!(res, u, order)
    length(u) == order + 1 || error("Length of u must be equal to order + 1")

    if order == 1
        # {0.5, 0.5}
        @. res = (u[1] + u[2]) * 0.5
    elseif order == 2
        # {-0.125, 0.75, 0.375}
        @. res = -0.125 * u[1] + 0.75 * u[2] + 0.375 * u[3]
    elseif order == 3
        # {0.0625, -0.3125, 0.9375, 0.3125}
        @. res = 0.0625 * u[1] - 0.3125 * u[2] + 0.9375 * u[3] + 0.3125 * u[4]
    elseif order == 4
        # {-0.0390625, 0.21875, -0.546875, 1.09375, 0.2734375}
        @. res =
            -0.0390625 * u[1] + 0.21875 * u[2] - 0.546875 * u[3] +
            1.09375 * u[4] +
            0.2734375 * u[5]
    else
        println("Time interpolation order not supported yet: order = ", order)
        exit()
    end
end

#===============================================================================
Functions needed by Mongwane's subcycling method
===============================================================================#
function calc_kfs_from_kcs!(kfs, kcs, dtc, tmp, interp_in_time::Bool)
    theta = interp_in_time ? 0.5 : 0.0
    dtf = 0.5 * dtc
    d1yc = tmp[1]
    d2yc = tmp[2]
    d3yc = tmp[3]
    rk4_dense_output_dy!(d1yc, theta, dtc, kcs)
    rk4_dense_output_d2y!(d2yc, theta, dtc, kcs)
    rk4_dense_output_d3y!(d3yc, theta, dtc, kcs)

    fyd2yc = 4 * (kcs[3] - kcs[2]) / dtc^3

    @. kfs[1] = dtf * d1yc
    @. kfs[2] = dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc - fyd2yc)
    @. kfs[3] = dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc + fyd2yc)

    return nothing
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
                fidx = if (j == 1)
                    i + num_buffer_points
                else
                    num_total_points - i + 1 - num_buffer_points
                end
                cidx = parent_map[fidx]
                w = transition_profile(a, b, grid.levels[l].x[fidx])
                if is_aligned[fidx]
                    kcs = [coarse_level.k[m][v][cidx] for m in 1:4]
                    ys = interp_in_time ? DenseOutput.y(0.5, uc_p[cidx], kcs) : uc_p[cidx]
                    uf[fidx] = (1 - w) * ys + w * uf[fidx]
                else
                    nys = spatial_interpolation_order + 1
                    soffset = (mod(nys, 2) == 0) ? div(nys, 2) : div(nys, 2) + 1
                    ys = zeros(Float64, nys)
                    for ic in 1:nys
                        ic_grid = cidx + ic - soffset
                        kcs = [coarse_level.k[m][v][ic_grid] for m in 1:4]
                        ys[ic] = if interp_in_time
                            DenseOutput.y(0.5, uc_p[ic_grid], kcs)
                        else
                            uc_p[ic_grid]
                        end
                    end
                    uf[fidx] =
                        (1 - w) * interpolate(ys, soffset, spatial_interpolation_order) +
                        w * uf[fidx]
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
function prolongation_mongwane!(grid, l)
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_total_points,
        num_buffer_points,
        spatial_interpolation_order,
        time_interpolation_order,
    ) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    interp_in_time = time_interpolation_order > 0

    dtc = coarse_level.dt

    statef = fine_level.state
    statec = coarse_level.state

    kf = fine_level.k
    kc = coarse_level.k

    uf = statef[end]
    uc_p = statec[end - 1]

    # j: 1: left, 2: right
    for i in 1:num_buffer_points, j in 1:2
        fidx = if (j == 1)
            (num_buffer_points + 1 - i)
        else
            (num_total_points - num_buffer_points + i)
        end
        is_aligned = mod(fidx, 2) == 0

        if is_aligned
            cidx = fidx2cidx(fine_level, fidx)

            kcs = [view(kc[m], cidx, :) for m in 1:4]
            kfs = [view(kf[m], fidx, :) for m in 1:3]
            calc_kfs_from_kcs!(kfs, kcs, dtc, interp_in_time)
            # setting u
            uf[fidx] = interp_in_time ? DenseOutput.y(0.5, uc_p[cidx], kcs) : uc_p[cidx]
        else
            cidx = fidx2cidx(fine_level, fidx - 1)

            kfss = zeros(Float64, 3, num_spatial_interpolation_points)
            ys = zeros(Float64, num_spatial_interpolation_points)
            for ic in 1:nys
                ic_grid = cidx + ic - soffset
                kcs = [coarse_level.k[m][v][ic_grid] for m in 1:4]
                kfss[:, ic] = calc_kfs_from_kcs!(kcs, dtc, interp_in_time)
                ys[ic] = if interp_in_time
                    DenseOutput.y(0.5, uc_p[ic_grid], kcs)
                else
                    uc_p[ic_grid]
                end
            end
            # setting k
            for m in 1:3
                fine_level.k[m][v][fidx] = interpolate(
                    kfss[m, :], soffset, spatial_interpolation_order
                )
            end
            # setting u
            uf[fidx] = interpolate(ys, soffset, spatial_interpolation_order)
        end
    end
end

#===============================================================================
prolongation!:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
    * time interpolation is used when time_interpolation_order > 0
===============================================================================#
function prolongation!(
    grid::Grid{NumState,NumDiagnostic}, l::Int
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_total_points,
        num_buffer_points,
        spatial_interpolation_order,
        time_interpolation_order,
    ) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    interp_in_time = time_interpolation_order > 0

    statef = fine_level.state
    statec = coarse_level.state

    uf = statef[end]
    uc_p = statec[end - 1]

    # j: 1: left, 2: right
    if interp_in_time
        buffer = zeros(Float64, num_spatial_interpolation_points, NumState)
        for i in 1:num_buffer_points, j in 1:2
            # from nearest point to the boundary
            fidx = if (j == 1)
                (num_buffer_points + 1 - i)
            else
                (num_total_points - num_buffer_points + i)
            end
            is_aligned = mod(fidx, 2) == 0
            if is_aligned
                cidx = fidx2cidx(fine_level, fidx)
                # time interpolation
                prolongation_time_interpolate!(
                    @view(uf[fidx, :]),
                    [@view(statec[m][cidx, :]) for m in 1:(time_interpolation_order + 1)],
                    time_interpolation_order,
                )
            else
                cidx = fidx2cidx(fine_level, fidx - 1)
                for ic in 1:num_spatial_interpolation_points
                    ic_grid = cidx + ic - soffset
                    # time interpolation
                    prolongation_time_interpolate!(
                        @view(buffer[ic, :]),
                        [
                            @view(statec[m][ic_grid, :]) for
                            m in 1:(time_interpolation_order + 1)
                        ],
                        time_interpolation_order,
                    )
                end
                # spatial interpolation
                prolongation_spatial_interpolate!(
                    @view(uf[fidx, :]),
                    [@view(buffer[m, :]) for m in 1:num_spatial_interpolation_points],
                    soffset,
                    spatial_interpolation_order,
                )
            end
        end
    else
        for i in 1:num_buffer_points, j in 1:2
            fidx = if (j == 1)
                (num_buffer_points + 1 - i)
            else
                (num_total_points - num_buffer_points + i)
            end
            is_aligned = mod(fidx, 2) == 0
            if is_aligned
                cidx = fidx2cidx(fine_level, fidx)
                uf[fidx] = uc_p[cidx]
            else
                cidx = fidx2cidx(fine_level, fidx - 1)
                # spatial interpolation
                prolongation_spatial_interpolate!(
                    @view(uf[fidx, :]),
                    [
                        @view(uc_p[cidx + ic - soffset, :]) for
                        ic in 1:num_spatial_interpolation_points
                    ],
                    soffset,
                    spatial_interpolation_order,
                )
            end
        end
    end
end
