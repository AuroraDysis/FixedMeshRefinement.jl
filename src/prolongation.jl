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
        error("Spatial interpolation order not supported yet: order = ", order)
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
        error("Time interpolation order not supported yet: order = ", order)
    end
end

#===============================================================================
Functions needed by Mongwane's subcycling method
===============================================================================#
function calc_Yn_from_kcs!(Yn_buffer, yn, kcs, dtc, interp_in_time::Bool, dytmp)
    theta = interp_in_time ? 0.5 : 0.0

    d1y = dytmp[1]
    d2y = dytmp[2]
    d3y = dytmp[3]

    Y1 = Yn_buffer[1]
    Y2 = Yn_buffer[2]
    Y3 = Yn_buffer[3]
    Y4 = Yn_buffer[4]

    if interp_in_time
        rk4_dense_output_y!(Y1, theta, dtc, yn, kcs)
    else
        Y1 .= yn
    end

    rk4_dense_output_dy!(d1y, theta, dtc, yn, kcs)
    rk4_dense_output_d2y!(d2y, theta, dtc, yn, kcs)
    rk4_dense_output_d3y!(d3y, theta, dtc, yn, kcs)

    # ratio between coarse and fine cell size (2 to 1 MR case)

    r = 0.5
    r2 = r * r
    k2 = kcs[2]
    k3 = kcs[3]
    dtf = 0.5 * dtc
    @. Y2 = Y1 + dtf * r * d1y
    @. Y3 = Y1 + dtf * (0.5 * d1y + 0.25 * r * d2y + 0.0625 * r2 * (d3y + 4.0 * (k3 - k2)))
    @. Y4 = Y1 + dtf * (r * d1y + 0.5 * r * d2y + 0.125 * r2 * (d3y - 4.0 * (k3 - k2)))

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
        error("Transition profile (type $type) is not supported yet.")
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

    kc = coarse_level.k

    uf = statef[end]
    uc_p = statec[end - 1]

    aligned_buffer = zeros(Float64, NumState)
    spatial_buffer = zeros(Float64, num_spatial_interpolation_points, NumState)

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
                kcs = [view(kc[m], cidx, :) for m in 1:4]
                yn = @view(uc_p[cidx, :])
                if interp_in_time
                    rk4_dense_output_y!(aligned_buffer, 0.5, dtc, yn, kcs)
                else
                    aligned_buffer .= yn
                end
            else
                cidx = fidx2cidx(fine_level, fidx - 1)
                for ic in 1:num_spatial_interpolation_points
                    ic_grid = cidx + ic - soffset
                    kcs = [view(kc[m], ic_grid, :) for m in 1:4]
                    yn = @view(uc_p[ic_grid, :])
                    if interp_in_time
                        rk4_dense_output_y!(@view(spatial_buffer[ic, :]), 0.5, dtc, yn, kcs)
                    else
                        spatial_buffer[ic, :] .= yn
                    end
                end
                prolongation_spatial_interpolate!(
                    aligned_buffer,
                    [
                        @view(spatial_buffer[ic, :]) for
                        ic in 1:num_spatial_interpolation_points
                    ],
                    soffset,
                    spatial_interpolation_order,
                )
            end
            @. uf[fidx, :] = (1 - w) * aligned_buffer + w * @view(uf[fidx, :])
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

    # Yn buffer
    dytmp = [zeros(Float64, NumState) for _ in 1:3]

    # dir: 1: left, 2: right
    for dir in 1:2, i in 1:num_buffer_points
        fidx = buffer_indices[dir][i]
        is_aligned = mod(fidx, 2) == 0

        if is_aligned
            cidx = fidx2cidx(fine_level, fidx)
            buffer = [view(Yn_buffer[rk_stage], i, :, dir) for rk_stage in 1:4]
            yn = @view(uc_p[cidx, :])
            kcs = [view(kc[m], cidx, :) for m in 1:4]
            calc_Yn_from_kcs!(buffer, yn, kcs, dtc, interp_in_time, dytmp)
        else
            cidx = fidx2cidx(fine_level, fidx - 1)

            for ic in 1:num_spatial_interpolation_points
                ic_grid = cidx + ic - soffset
                kcs = [view(kc[m], ic_grid, :) for m in 1:4]
                yn = @view(uc_p[ic_grid, :])
                sbuffer = [view(spatial_buffer[rk_stage], ic, :) for rk_stage in 1:4]
                calc_Yn_from_kcs!(sbuffer, yn, kcs, dtc, interp_in_time, dytmp)
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
