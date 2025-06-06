export Level, Grid

mutable struct Level{NumState,NumDiagnostic}
    num_interior_points::Int64  # num of interior grid points
    num_ghost_points::Int64  # num of ghost points
    num_buffer_points::Int64  # num of buffer points
    num_transition_points::Int64  # num of transition zone points
    num_total_points::Int64  # num of all grid points
    finite_difference_order::Int64  # finite difference order
    spatial_interpolation_order::Int64  # interpolation order in space
    domain_box::Vector{Float64}  # size computational domain (interior)
    spatial_step::Float64
    time_step::Float64
    time::Float64
    dissipation::Float64
    is_base_level::Bool
    parent_map::Vector{Int64}  # map between indexes of current and its parent level
    is_aligned::Vector{Bool}  # if grid aligned with coarse grid

    # data
    coordinates::Vector{Float64}  # coordinates
    state::Vector{Vector{Float64}}  # state vectors
    state_prev::Vector{Vector{Float64}}  # previous state vectors
    state_prev_prev::Vector{Vector{Float64}}  # previous previous state vectors
    rhs::Vector{Vector{Float64}}  # rhs of state vectors
    workspace::Vector{Vector{Float64}}  # intermediate state vectors
    # intermediate state vectors for new subcycling
    runge_kutta_stages::Vector{Vector{Vector{Float64}}}

    # diagnostic variables
    diag_state::Vector{Vector{Float64}}  # state vectors for diagnostic variables

    function Level{NumState,NumDiagnostic}(
        num_interior_points,
        num_ghost_points,
        num_buffer_points,
        num_transition_points,
        finite_difference_order,
        spatial_interpolation_order,
        domain_box,
        dt,
        t,
        dissipation,
        is_base_level,
        parent_map,
        is_aligned,
    )
        num_total_points = num_interior_points + 2 * num_buffer_points
        spatial_step = (domain_box[2] - domain_box[1]) / (num_interior_points - 1)

        noffset = (num_total_points - num_interior_points) / 2  # take account of buffer zone
        xmin = domain_box[1] - noffset * spatial_step
        xmax = domain_box[2] + noffset * spatial_step
        coordinates = collect(LinRange(xmin, xmax, num_total_points))
        state = Vector{Vector{Float64}}(undef, NumState)
        state_prev = Vector{Vector{Float64}}(undef, NumState)
        state_prev_prev = Vector{Vector{Float64}}(undef, NumState)
        rhs = Vector{Vector{Float64}}(undef, NumState)
        workspace = Vector{Vector{Float64}}(undef, NumState)
        runge_kutta_stages = Vector{Vector{Vector{Float64}}}(undef, 4)
        diag_state = Vector{Vector{Float64}}(undef, NumDiagnostic)
        for j in 1:4
            runge_kutta_stages[j] = Vector{Vector{Float64}}(undef, NumState)
        end

        for i in 1:NumState
            state[i] = fill(NaN, num_total_points)
            state_prev[i] = fill(NaN, num_total_points)
            state_prev_prev[i] = fill(NaN, num_total_points)
            rhs[i] = fill(NaN, num_total_points)
            workspace[i] = fill(NaN, num_total_points)
            for j in 1:4
                runge_kutta_stages[j][i] = fill(NaN, num_total_points)
            end
        end

        for i in 1:NumDiagnostic
            diag_state[i] = fill(NaN, num_total_points)
        end

        return new(
            num_interior_points,
            num_ghost_points,
            num_buffer_points,
            num_transition_points,
            num_total_points,
            finite_difference_order,
            spatial_interpolation_order,
            domain_box,
            spatial_step,
            dt,
            t,
            dissipation,
            is_base_level,
            parent_map,
            is_aligned,
            # data
            coordinates,
            state,
            state_prev,
            state_prev_prev,
            rhs,
            workspace,
            runge_kutta_stages,
            diag_state,
        )
    end
end

mutable struct Grid{NumState,NumDiagnostic}
    levels::Vector{Level{NumState,NumDiagnostic}}
    base_dt::Float64
    time::Float64
    use_subcycling::Bool  # turn on subcycling or not

    function Grid{NumState,NumDiagnostic}(
        base_level_num_points,  # num of interior grid points at base level
        domain_boxes::Vector{Vector{Float64}},
        num_ghost_points,
        num_buffer_points;
        num_transition_points=3,
        finite_difference_order=4,
        spatial_interpolation_order=5,
        cfl_number=0.25,
        initial_time=0.0,
        dissipation=0.0,
        use_subcycling=true,
    )
        # build the first level (base level)
        base_level_spatial_step =
            (domain_boxes[1][2] - domain_boxes[1][1]) / (base_level_num_points - 1)
        base_dt = if use_subcycling
            cfl_number * base_level_spatial_step
        else
            cfl_number * base_level_spatial_step / 2^(length(domain_boxes) - 1)
        end
        base_level = Level{NumState,NumDiagnostic}(
            base_level_num_points,
            num_ghost_points,
            num_buffer_points,
            num_transition_points,
            finite_difference_order,
            spatial_interpolation_order,
            domain_boxes[1],
            base_dt,
            initial_time,
            dissipation,
            true,
            [],
            [],
        )
        levels = Vector{Level{NumState,NumDiagnostic}}([base_level])
        # build the rest levels
        for i in 2:length(domain_boxes)
            level_spatial_step = base_level_spatial_step / 2^(i - 1)
            level_time_step = (use_subcycling ? cfl_number * level_spatial_step : base_dt)
            parent_level = levels[i - 1]  # level lower than the current level (parent level)
            # if we refine parent level everywhere
            parent_level_grid_points = LinRange(
                parent_level.domain_box[1],
                parent_level.domain_box[2],
                (parent_level.num_interior_points - 1) * 2 + 1,
            )
            # find those two which are closest to current level boundaries
            min_index = argmin(abs.(parent_level_grid_points .- domain_boxes[i][1]))
            max_index = argmin(abs.(parent_level_grid_points .- domain_boxes[i][2]))
            #imin = findall(x->abs(x - domain_boxes[i][1]) <= level_spatial_step + 1e-12, parent_level_grid_points)[1]
            #imax = findall(x->abs(x - domain_boxes[i][2]) <= level_spatial_step + 1e-12, parent_level_grid_points)[end]
            level_domain = [
                parent_level_grid_points[min_index], parent_level_grid_points[max_index]
            ]
            level_num_interior_points = (max_index - min_index) + 1  # (floor(Int, (parent_level_grid_points[max_index] - parent_level_grid_points[min_index]) / level_spatial_step)) + 1
            # maps between two levels
            parent_map =
                div.(
                    (
                        ((min_index - num_buffer_points):(max_index + num_buffer_points)) .+
                        1
                    ),
                    2,
                ) .+ num_buffer_points
            is_aligned =
                mod.(
                    (
                        ((min_index - num_buffer_points):(max_index + num_buffer_points)) .+
                        1
                    ),
                    2,
                ) .== 0
            # build level
            push!(
                levels,
                Level{NumState,NumDiagnostic}(
                    level_num_interior_points,
                    num_ghost_points,
                    num_buffer_points,
                    num_transition_points,
                    finite_difference_order,
                    spatial_interpolation_order,
                    level_domain,
                    level_time_step,
                    initial_time,
                    dissipation,
                    false,
                    parent_map,
                    is_aligned,
                ),
            )
        end

        # construct
        return new{NumState,NumDiagnostic}(levels, base_dt, initial_time, use_subcycling)
    end
end

function Base.show(
    io::IO, grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    println(io, "Grid Structure:")
    println(io, "  use_subcycling = ", grid.use_subcycling)
    for i in 1:length(grid.levels)
        println(io, "level[", i, "],")
        println(io, "  num_interior_points       = ", grid.levels[i].num_interior_points)
        println(io, "  num_ghost_points          = ", grid.levels[i].num_ghost_points)
        println(io, "  num_buffer_points         = ", grid.levels[i].num_buffer_points)
        println(io, "  num_transition_points     = ", grid.levels[i].num_transition_points)
        println(
            io, "  finite_difference_order   = ", grid.levels[i].finite_difference_order
        )
        println(
            io,
            "  spatial_interpolation_order = ",
            grid.levels[i].spatial_interpolation_order,
        )
        if length(grid.levels[i].parent_map) == grid.levels[i].num_total_points
            println(
                io,
                "  ibox   = ",
                [
                    grid.levels[i].parent_map[1 + grid.levels[i].num_buffer_points],
                    grid.levels[i].parent_map[grid.levels[i].num_total_points - grid.levels[i].num_buffer_points],
                ],
                ", ",
                [
                    grid.levels[i].is_aligned[1 + grid.levels[i].num_buffer_points],
                    grid.levels[i].is_aligned[grid.levels[i].num_total_points - grid.levels[i].num_buffer_points],
                ],
            )
        end
        println(io, "  domain_box    = ", grid.levels[i].domain_box)
        println(io, "  spatial_step  = ", grid.levels[i].spatial_step)
        println(io, "  time_step     = ", grid.levels[i].time_step)
        println(io, "  time          = ", grid.levels[i].time)
        println(io, "  dissipation   = ", grid.levels[i].dissipation)
    end
end
