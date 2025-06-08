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
function calc_Yn_from_kcs!(Yn_buffer, yn, kcs, dtc, interp_in_time::Bool)
    theta = interp_in_time ? 0.5 : 0.0
    dtf = 0.5 * dtc

    Y1 = Yn_buffer[1]
    Y2 = Yn_buffer[2]
    Y3 = Yn_buffer[3]
    Y4 = Yn_buffer[4]

    if interp_in_time
        rk4_dense_output_y!(@view(Y1[i, :, dir]), 0.5, dtc, yn, kcs)
    else
        Y1 .= yn
    end

    rk4_dense_output_dy!(Y2, theta, dtc, yn, kcs)
    rk4_dense_output_d2y!(Y3, theta, dtc, yn, kcs)
    rk4_dense_output_d3y!(Y4, theta, dtc, yn, kcs)

    fyd2yc = 4 * (kcs[3] - kcs[2]) / dtc^3

    @. Y2[1] = yn + dtf * d1yc
    @. Y3[2] = dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc - fyd2yc)
    @. Y4[3] = dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc + fyd2yc)

    return nothing
end

function transition_profile(xl, xh, x; type=1)
    t0 = (x - xl) / (xh - xl)
    t = if t0 < 0.0
        0.0
    elseif t0 > 1.0
        1.0
    else
        t0
    end

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
    ) = fine_level

    # for transition zone
    domain_box = fine_level.domain_box
    dxf = fine_level.dx

    isapprox(domain_box[1], fine_level.x[1 + num_buffer_points]; rtol=1e-12) ||
        error("domain_box[1] != fine_level.x[1 + num_buffer_points]")
    isapprox(
        domain_box[2], fine_level.x[num_total_points - num_buffer_points]; rtol=1e-12
    ) || error("domain_box[2] != fine_level.x[num_total_points - num_buffer_points]")

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    statef = fine_level.state
    statec = coarse_level.state

    uf = statef[end]
    uc_p = statec[end - 1]

    for dir in 1:2  # left or right
        a = dir == 1 ? domain_box[1] : domain_box[2]
        b = if dir == 1
            domain_box[1] + (num_transition_points - 1) * dxf
        else
            domain_box[2] - (num_transition_points - 1) * dxf
        end

        for i in 1:num_transition_points
            fidx = if dir == 1
                num_buffer_points + i
            else
                num_total_points - num_buffer_points + 1 - i
            end
            is_aligned = mod(i + 1, 2) == 0
            w = transition_profile(a, b, fine_level.x[fidx])
            if is_aligned
                cidx = fidx2cidx(fine_level, fidx)
                kcs = [coarse_level.k[m][v][cidx] for m in 1:4]
                ys = interp_in_time ? DenseOutput.y(0.5, uc_p[cidx], kcs) : uc_p[cidx]
                uf[fidx] = (1 - w) * ys + w * uf[fidx]
            else
                cidx = fidx2cidx(fine_level, fidx - 1)
                ys = zeros(Float64, num_spatial_interpolation_points)
                for ic in 1:num_spatial_interpolation_points
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

#===============================================================================
prolongation_mongwane!: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function prolongation_mongwane!(grid, l, interp_in_time::Bool)
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (; num_buffer_points, spatial_interpolation_order, buffer_indices) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    dtc = coarse_level.dt

    statef = fine_level.state
    statec = coarse_level.state

    kc = coarse_level.k

    uf = statef[end]
    uc_p = statec[end - 1]

    # buffer points for Yn
    Yn_buffer = fine_level.Yn_buffer

    # buffer points for spatial interpolation
    spatial_buffer = [
        zeros(Float64, num_spatial_interpolation_points, NumState) for _ in 1:4
    ]

    # dir: 1: left, 2: right
    for dir in 1:2, i in 1:num_buffer_points
        fidx = buffer_indices[dir][i]
        is_aligned = mod(fidx, 2) == 0

        if is_aligned
            cidx = fidx2cidx(fine_level, fidx)
            buffer = [view(Yn_buffer[rk_stage], i, :, dir) for rk_stage in 1:4]
            yn = @view(uc_p[cidx, :])
            kcs = [view(kc[m], cidx, :) for m in 1:4]
            calc_Yn_from_kcs!(buffer, yn, kcs, dtc, interp_in_time)
        else
            cidx = fidx2cidx(fine_level, fidx - 1)

            for ic in 1:num_spatial_interpolation_points
                ic_grid = cidx + ic - soffset
                kcs = [view(kc[m], ic_grid, :) for m in 1:4]
                yn = @view(uc_p[ic_grid, :])
                sbuffer = [view(spatial_buffer[rk_stage], ic, :) for rk_stage in 1:4]
                calc_Yn_from_kcs!(sbuffer, yn, kcs, dtc, interp_in_time)
            end
            for rk_stage in 1:4
                prolongation_spatial_interpolate!(
                    view(Yn_buffer[rk_stage], i, :, dir),
                    [
                        view(spatial_buffer[rk_stage], ic, :) for
                        ic in 1:num_spatial_interpolation_points
                    ],
                    soffset,
                    spatial_interpolation_order,
                )
            end
        end
    end

    # fill fine level buffer
    fill_buffer!(uf, fine_level, 1)

    return nothing
end

#===============================================================================
prolongation!:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
    * time interpolation is used when time_interpolation_order > 0
===============================================================================#
function prolongation!(
    grid::Grid{NumState,NumDiagnostic}, l::Int, interp_in_time::Bool
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_buffer_points,
        spatial_interpolation_order,
        time_interpolation_order,
        buffer_indices,
    ) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    statef = fine_level.state
    statec = coarse_level.state

    uf = statef[end]
    uc_p = statec[end - 1]

    if interp_in_time
        buffer = zeros(Float64, num_spatial_interpolation_points, NumState)
        # dir: 1: left, 2: right
        for dir in 1:2, i in 1:num_buffer_points
            fidx = buffer_indices[dir][i]
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
        # dir: 1: left, 2: right
        for dir in 1:2, i in 1:num_buffer_points
            fidx = buffer_indices[dir][i]
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
