export prolongate!, prolongate_mongwane!, apply_transition_zone!

"""
    spatial_interpolate!(res, u, i, order)

Performs spatial interpolation for prolongation using a set of predefined stencils of a given `order`.
The result is stored in `res`. The interpolation uses values from `u` around index `i`.

# Arguments
- `res`: The destination array for the interpolated values.
- `u`: The source data array.
- `i::Int`: The index in `u` to center the interpolation stencil.
- `order::Int`: The order of spatial interpolation (1 through 5 are supported).
"""
function spatial_interpolate!(res, u, i, order)
    if order == 1
        # {1/2, 1/2}
        @.. res = (u[i] + u[i + 1]) * (1 / 2)
    elseif order == 2
        # {-(1/8), 3/4, 3/8}
        @.. res = -(1 / 8) * u[i - 1] + (3 / 4) * u[i] + (3 / 8) * u[i + 1]
    elseif order == 3
        # {-(1/16),9/16,9/16,-(1/16)}
        @.. res = -(1 / 16) * (u[i - 1] + u[i + 2]) + (9 / 16) * (u[i] + u[i + 1])
    elseif order == 4
        # {3/128, -(5/32), 45/64, 15/32, -(5/128)}
        @.. res =
            (3 / 128) * u[i - 2] - (5 / 32) * u[i - 1] +
            (45 / 64) * u[i] +
            (15 / 32) * u[i + 1] - (5 / 128) * u[i + 2]
    elseif order == 5
        # {3/256, -25/256, 75/128, 75/128, -25/256, 3/256}
        @.. res =
            (3 / 256) * (u[i - 2] + u[i + 3]) - (25 / 256) * (u[i - 1] + u[i + 2]) +
            (75 / 128) * (u[i] + u[i + 1])
    else
        error("Spatial interpolation order not supported yet: order = ", order)
    end
end

"""
    time_interpolate!(res, u, order)

Performs time interpolation for prolongation. `u` should be ordered from oldest to newest.
The result is stored in `res`.

# Arguments
- `res`: The destination array for the interpolated values.
- `u`: A vector of states at different time levels, ordered from oldest to newest.
- `order::Int`: The order of time interpolation (1 through 4 are supported).
"""
# u from oldest to newest
function time_interpolate!(res, u, order)
    length(u) == order + 1 || error("Length of u must be equal to order + 1")

    if order == 1
        # {1/2, 1/2}
        @.. res = (u[1] + u[2]) * (1 / 2)
    elseif order == 2
        # {-(1/8), 3/4, 3/8}
        @.. res = -(1 / 8) * u[1] + (3 / 4) * u[2] + (3 / 8) * u[3]
    else
        error("Time interpolation order not supported yet: order = ", order)
    end
end

#===============================================================================
Functions needed by Mongwane's subcycling method
===============================================================================#
"""
    rk4_dense_output_y!(y, theta, h, yn, k)

Calculate the dense output for a 4th-order Runge-Kutta method at a fractional time
`theta` within a time step `h`. This is used for interpolation in time within a single
time step.

# Arguments
- `y`: The output array for the interpolated state.
- `theta::Float64`: The fractional time within the interval `[0, 1]`.
- `h::Float64`: The time step size.
- `yn`: The state at the beginning of the time step.
- `k`: A vector of RK4 stages `k1` through `k4`.
"""
function rk4_dense_output_y!(y, theta, h, yn, k)
    theta2 = theta * theta
    theta3 = theta2 * theta
    ck1 = theta - 3 * theta2 / 2 + 2 * theta3 / 3
    ck23 = theta2 - 2 * theta3 / 3
    ck4 = -theta2 / 2 + 2 * theta3 / 3
    @.. y = yn + (ck1 * k[1] + ck23 * (k[2] + k[3]) + ck4 * k[4]) * h
    return nothing
end

"""
    calc_Yn_from_kcs!(Yn_buffer, yn, kcs, dtc, interp_in_time, dytmp)

Calculate the intermediate `Y` values needed for Mongwane's subcycling method. These
values represent the state on the coarse grid at intermediate times required by the
fine grid's time stepping.

# Arguments
- `Yn_buffer`: Buffer to store the calculated `Y` values.
- `yn`: State on the coarse grid at the beginning of the step.
- `kcs`: The RK4 stages from the coarse grid time step.
- `dtc::Float64`: The coarse grid time step.
- `interp_in_time::Bool`: Flag indicating if interpolation in time is needed.
- `dytmp`: Temporary storage for derivatives.
"""
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
        @.. Y1 = yn
    end

    theta2 = theta * theta
    let ck1 = 1 - 3 * theta + 2 * theta2,
        ck23 = 2 * theta - 2 * theta2,
        ck4 = -theta + 2 * theta2

        @.. d1y = ck1 * kcs[1] + ck23 * (kcs[2] + kcs[3]) + ck4 * kcs[4]
    end
    let ck1 = -3 + 4 * theta, ck23 = 2 - 4 * theta, ck4 = -1 + 4 * theta
        @.. d2y = (ck1 * kcs[1] + ck23 * (kcs[2] + kcs[3]) + ck4 * kcs[4]) * dtc
    end
    @.. d3y = 4 * (kcs[1] - kcs[2] - kcs[3] + kcs[4])

    # ratio between coarse and fine cell size (2 to 1 MR case)

    r = 0.5
    r2 = r * r
    k2 = kcs[2]
    k3 = kcs[3]
    dtf = r * dtc
    @.. Y2 = Y1 + 0.5 * dtf * d1y
    @.. Y3 = Y1 + dtf * (0.5 * d1y + 0.25 * r * d2y + 0.0625 * r2 * (d3y + 4.0 * (k3 - k2)))
    @.. Y4 = Y1 + dtf * (d1y + 0.5 * r * d2y + 0.125 * r2 * (d3y - 4.0 * (k3 - k2)))

    return nothing
end

"""
    transition_profile(xl, xh, x; type=1)

Compute a transition profile value at `x` between `xl` and `xh`. This is used to
smoothly blend solutions in transition zones.

# Arguments
- `xl::Float64`: The lower bound of the transition interval.
- `xh::Float64`: The upper bound of the transition interval.
- `x::Float64`: The point at which to evaluate the profile.
- `type::Int`: The type of transition profile (1: boxstep, 2: smoothstep, 3: smootherstep).
"""
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

"""
    apply_transition_zone!(grid::Grid, l::Int, interp_in_time::Bool)

Apply a transition zone to smoothly connect coarse and fine grid solutions at the
boundaries of a refinement level. This function blends the fine grid solution with
a prolonged solution from the coarse grid.

# Arguments
- `grid::Grid`: The grid structure.
- `l::Int`: The fine level index.
- `interp_in_time::Bool`: Whether to use time interpolation.
"""
function apply_transition_zone!(
    grid::Grid{NumState,NumDiagnostic}, l::Int, interp_in_time::Bool
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_interior_points,
        num_transition_points,
        is_physical_boundary,
        spatial_interpolation_order,
    ) = fine_level

    # for transition zone
    domain_box = fine_level.domain_box
    dxf = fine_level.dx

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    uf = get_state(fine_level)
    uc_p = get_state(coarse_level, -1)

    dtc = coarse_level.dt

    kc = [get_rk_stage(coarse_level, i) for i in 1:4]

    aligned_buffer = MVector{NumState,Float64}(undef)
    spatial_buffer = zeros(Float64, num_spatial_interpolation_points, NumState)

    fine_x = get_x(fine_level)

    # left or right
    for dir in 1:2
        # skip physical boundary
        if is_physical_boundary[dir]
            continue
        end

        a = domain_box[dir]
        b = if dir == 1
            domain_box[dir] + (num_transition_points[dir] - 1) * dxf
        else
            domain_box[dir] - (num_transition_points[dir] - 1) * dxf
        end

        for i in 1:num_transition_points[dir]
            fidx = if dir == 1
                i
            else
                num_interior_points + 1 - i
            end
            is_aligned = isodd(i)

            w = transition_profile(a, b, fine_x[fidx])
            if is_aligned
                cidx = fine_to_coarse_index(fine_level, fidx)
                kcs = [view(kc[m], cidx, :) for m in 1:4]
                yn = @view(uc_p[cidx, :])
                if interp_in_time
                    rk4_dense_output_y!(aligned_buffer, 0.5, dtc, yn, kcs)
                else
                    aligned_buffer .= yn
                end
            else
                cidx = fine_to_coarse_index(fine_level, fidx - 1)
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
                spatial_interpolate!(
                    aligned_buffer,
                    [
                        @view(spatial_buffer[ic, :]) for
                        ic in 1:num_spatial_interpolation_points
                    ],
                    soffset,
                    spatial_interpolation_order,
                )
            end
            @.. uf[fidx, :] = (1 - w) * aligned_buffer + w * @view(uf[fidx, :])
        end
    end
end

#===============================================================================
prolongate_mongwane!: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
"""
    prolongate_mongwane!(grid::Grid, l::Int, interp_in_time::Bool)

Perform prolongation from a coarse grid (`l-1`) to a fine grid (`l`) using Mongwane's
subcycling-in-time method. This is used to set the ghost cell data for the fine grid
when subcycling is enabled.

# Arguments
- `grid::Grid`: The grid structure.
- `l::Int`: The fine level index.
- `interp_in_time::Bool`: Whether to use time interpolation for the coarse grid state.
"""
function prolongate_mongwane!(
    grid::Grid{NumState,NumDiagnostic}, l::Int, interp_in_time::Bool
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (; num_boundary_points, spatial_interpolation_order, is_physical_boundary) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    dtc = coarse_level.dt

    kc = [get_rk_stage(coarse_level, i) for i in 1:4]

    uf = get_state(fine_level)
    uc_p = get_state(coarse_level, -1)

    boundary_indices = get_boundary_indices(fine_level)

    # buffer points for Yn
    Yn_buffer = fine_level.Yn_buffer

    # buffer points for spatial interpolation
    spatial_buffer = [
        zeros(Float64, num_spatial_interpolation_points, NumState) for _ in 1:4
    ]

    # Yn buffer
    dytmp = [MVector{NumState,Float64}(undef) for _ in 1:3]

    # dir: 1: left, 2: right
    for dir in 1:2
        if is_physical_boundary[dir]
            continue
        end

        for i in 1:num_boundary_points[dir]
            fidx = boundary_indices[dir][i]
            is_aligned = iseven(i)

            if is_aligned
                cidx = fine_to_coarse_index(fine_level, fidx)
                buffer = [view(Yn_buffer[rk_stage], i, :, dir) for rk_stage in 1:4]
                yn = @view(uc_p[cidx, :])
                kcs = [view(kc[m], cidx, :) for m in 1:4]
                calc_Yn_from_kcs!(buffer, yn, kcs, dtc, interp_in_time, dytmp)
            else
                cidx = fine_to_coarse_index(fine_level, fidx - 1)

                for ic in 1:num_spatial_interpolation_points
                    ic_grid = cidx + ic - soffset
                    kcs = [view(kc[m], ic_grid, :) for m in 1:4]
                    yn = @view(uc_p[ic_grid, :])
                    sbuffer = [view(spatial_buffer[rk_stage], ic, :) for rk_stage in 1:4]
                    calc_Yn_from_kcs!(sbuffer, yn, kcs, dtc, interp_in_time, dytmp)
                end
                for rk_stage in 1:4
                    spatial_interpolate!(
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
    end

    # fill fine level buffer
    fill_buffer!(uf, fine_level, 1)

    return nothing
end

#===============================================================================
prolongate!:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
    * time interpolation is used when time_interpolation_order > 0
===============================================================================#
"""
    prolongate!(grid::Grid, l::Int, interp_in_time::Bool)

Perform prolongation from a coarse grid (`l-1`) to a fine grid (`l`). This function
fills the ghost cells of the fine grid by interpolating the data from the coarse grid
in space and, optionally, in time.

# Arguments
- `grid::Grid`: The grid structure.
- `l::Int`: The fine level index.
- `interp_in_time::Bool`: Whether to use time interpolation.
"""
function prolongate!(
    grid::Grid{NumState,NumDiagnostic}, l::Int, interp_in_time::Bool
) where {NumState,NumDiagnostic}
    fine_level = grid.levels[l]
    coarse_level = grid.levels[l - 1]

    (;
        num_boundary_points,
        spatial_interpolation_order,
        time_interpolation_order,
        is_physical_boundary,
    ) = fine_level

    num_spatial_interpolation_points = spatial_interpolation_order + 1
    soffset = if mod(num_spatial_interpolation_points, 2) == 0
        div(num_spatial_interpolation_points, 2)
    else
        div(num_spatial_interpolation_points, 2) + 1
    end

    uf = get_state(fine_level)
    uc_p = get_state(coarse_level, -1)
    boundary_indices = get_boundary_indices(fine_level)

    buffer = zeros(Float64, num_spatial_interpolation_points, NumState)

    # dir: 1: left, 2: right
    for dir in 1:2
        if is_physical_boundary[dir]
            continue
        end

        for i in 1:num_boundary_points[dir]
            fidx = boundary_indices[dir][i]
            is_aligned = iseven(i)

            if interp_in_time
                if is_aligned
                    cidx = fine_to_coarse_index(fine_level, fidx)
                    time_data = [
                        view(get_state(coarse_level, m), cidx, :) for
                        m in (-time_interpolation_order):0
                    ]
                    # time interpolation
                    time_interpolate!(
                        @view(uf[fidx, :]), time_data, time_interpolation_order
                    )
                else
                    cidx = fine_to_coarse_index(fine_level, fidx - 1)
                    for ic in 1:num_spatial_interpolation_points
                        ic_grid = cidx + ic - soffset
                        time_data = [
                            view(get_state(coarse_level, m), ic_grid, :) for
                            m in (-time_interpolation_order):0
                        ]
                        # time interpolation
                        time_interpolate!(
                            @view(buffer[ic, :]), time_data, time_interpolation_order
                        )
                    end
                    # spatial interpolation
                    spatial_interpolate!(
                        @view(uf[fidx, :]),
                        [@view(buffer[m, :]) for m in 1:num_spatial_interpolation_points],
                        soffset,
                        spatial_interpolation_order,
                    )
                end
            else
                if is_aligned
                    cidx = fine_to_coarse_index(fine_level, fidx)
                    uf[fidx, :] .= @view(uc_p[cidx, :])
                else
                    cidx = fine_to_coarse_index(fine_level, fidx - 1)
                    # spatial interpolation
                    spatial_interpolate!(
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
end
