include("../../src/FixedMeshRefinement.jl")
using .FixedMeshRefinement

using Printf
using TOML
using FastBroadcast
using DelimitedFiles

include("boundary.jl")
include("initial_data.jl")
include("output_csv.jl")
include("output_hdf5.jl")
include("output_checkpoint.jl")
include("wave.jl")

function get(params, key, default)
    return haskey(params, key) ? params[key] : default
end

function nan_check(grid::Grid)
    has_nan = false
    num_levels = get_num_levels(grid)
    for l in 1:num_levels
        level = get_level(grid, l)
        u = get_state(level)
        (; num_boundary_points, parent_indices, is_physical_boundary) = level
        interior_indices = get_interior_indices(level)
        if any(isnan.(u[:, interior_indices]))
            has_nan = true
            # print all the nan indexes
            println("level $l:")
            println("  x: ", get_x(level))
            println("  num_boundary_points: ", num_boundary_points)
            println("  parent_indices: ", parent_indices)
            println("  is_physical_boundary: ", is_physical_boundary)
            for i in 1:(level.num_state_variables)
                println(
                    "  nan_points $(i): ", interior_indices[isnan.(u[i, interior_indices])]
                )
            end
        end
    end
    return has_nan
end

function main(params, out_dir; grid=nothing, start_step=0)
    ########################
    # Read Parameter Files #
    ########################
    stop_time::Float64 = get(params, "stop_time", -1.0)
    max_step::Int = get(params, "max_step", -1)
    out_every_dt_0d::Float64 = get(params, "out_every_0d", 0.01)
    out_every_dt_1d::Float64 = get(params, "out_every_1d", 0.1)
    cfl::Float64 = get(params, "cfl", 0.25)
    dissipation::Float64 = get(params, "dissipation", 0.0)
    subcycling::Bool = get(params, "subcycling", true)
    mongwane::Bool = get(params, "mongwane", false)
    apply_trans_zone::Bool = get(params, "apply_trans_zone", true)
    checkpoint_every_dt::Float64 = get(params, "checkpoint_every", 2.0)

    println("Parameters:")
    println("  cfl              = ", cfl)
    println("  mongwane         = ", mongwane)
    println("  trans_zone       = ", apply_trans_zone)
    println("  max_step         = ", max_step)
    println("  stop_time        = ", stop_time)
    println("  out_every_dt_0d  = ", out_every_dt_0d)
    println("  out_every_dt_1d  = ", out_every_dt_1d)
    println("  checkpoint_every_dt = ", checkpoint_every_dt)
    println("  out_dir          = ", out_dir)

    ########################
    # build grid structure #
    ########################
    p = (; dissipation, termination=Ref(false))
    if isnothing(grid)
        num_state_variables = 2
        num_diagnostic_variables = 1
        num_tmp_variables = 3

        num_interior_points = params["num_interior_points"]
        num_ghost_points = params["num_ghost_points"]
        num_buffer_points = get(params, "num_buffer_points", 4 * num_ghost_points)
        params_domain_boxes = params["domain_boxes"]
        domain_boxes = [NTuple{2,Float64}(box) for box in params_domain_boxes]
        num_transition_points = get(params, "num_transition_points", 3)
        spatial_interpolation_order = get(params, "spatial_interpolation_order", 5)
        initial_data = get(params, "initial_data", "gaussian")

        grid = Grid(
            num_state_variables,
            num_interior_points,
            domain_boxes,
            num_ghost_points,
            num_buffer_points;
            num_diagnostic_variables=num_diagnostic_variables,
            num_tmp_variables=num_tmp_variables,
            num_transition_points=num_transition_points,
            spatial_interpolation_order=spatial_interpolation_order,
            cfl=cfl,
            subcycling=subcycling,
        )

        # just for testing, if all levels are aligned with the physical boundary, then we excise some grid points
        # if all([level.is_physical_boundary[1] for level in grid.levels])
        #     shift_grid_boundaries!(grid, (-2, 0))
        # elseif all([level.is_physical_boundary[2] for level in grid.levels])
        #     shift_grid_boundaries!(grid, (0, -2))
        # end

        ###############
        # Intial Data #
        ###############
        println("Setting up initial conditions...")
        println("  initial data type: $initial_data")
        if initial_data == "gaussian"
            gaussian!(grid)
        elseif initial_data == "sinusoidal"
            sinusoidal!(grid)
        else
            error("Initial data type '$initial_data' unsupported yet")
        end

        apply_reflective_boundary_condition!(grid)
        if !mongwane
            march_backwards!(grid, p)
            apply_reflective_boundary_condition!(grid)
        end
    else
        println("Restarting from step $(start_step-1)")
    end

    # print grid
    num_levels = get_num_levels(grid)
    for l in 1:num_levels
        show(stdout, MIME("text/plain"), get_level(grid, l))
    end

    # find output directories
    i = 1
    begin
        out_dir_tmp = ""
        while i <= 100
            out_dir_tmp = joinpath(out_dir, "out_$(lpad(i, 2, '0'))")
            if !isdir(out_dir_tmp)
                break
            end
            i += 1
        end
        if i > 100
            error("Too many existing output directories, something might be wrong")
        end
        out_dir = out_dir_tmp
    end

    # create a folder for csv
    csv_dir = joinpath(out_dir, "csv")
    if !isdir(csv_dir)
        mkpath(csv_dir)
    end

    out_csv = OutputCSV(joinpath(csv_dir, "output.csv"), ["t", "Ebase", "E"])

    begin
        num_levels = get_num_levels(grid)

        for l in 1:num_levels
            level = get_level(grid, l)
            x = get_x(level)

            # Create filename with level number padded to 2 digits
            csv_filename = "x_$(lpad(l, 2, '0')).csv"
            csv_path = joinpath(csv_dir, csv_filename)
            writedlm(csv_path, x)
        end
    end

    data_dir = joinpath(out_dir, "data")
    if !isdir(data_dir)
        mkpath(data_dir)
    end
    out_h5 = OutputHDF5(data_dir, "data", grid)

    checkpoint_dir = joinpath(out_dir, "checkpoints")
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    Ebase, E = wave_energy(grid)
    E0 = E
    @printf(
        "t = %.4f, iteration %d. Ebase = %.17f, E = %.17f, diff = %.3g\n",
        grid.t,
        start_step - 1,
        Ebase,
        E,
        Ebase - E
    )

    if start_step == 1
        write_row(out_csv, (grid.t, Ebase, E))
        append_data(out_h5, grid, 1)
    end

    ##########
    # Evolve #
    ##########
    println("Start evolution...")

    runtime_start = time()
    slurm_job_end_time = if haskey(ENV, "SLURM_JOB_END_TIME")
        parse(Float64, ENV["SLURM_JOB_END_TIME"])
    else
        Inf
    end

    step = start_step
    out_every_0d = round(Int, out_every_dt_0d / grid.base_dt)
    out_every_1d = round(Int, out_every_dt_1d / grid.base_dt)
    start_t = grid.t
    next_checkpoint_t = if start_step == 0
        checkpoint_every_dt
    else
        start_t - mod(start_t, checkpoint_every_dt) + checkpoint_every_dt
    end

    while true
        if out_every_0d > 0 && mod(step, out_every0d) == 0
            Ebase, E = wave_energy(grid)
            @printf(
                "t = %.4f, iteration %d. dEbase = %.5g, dE = %.5g\n",
                grid.t,
                step,
                Ebase - E0,
                E - E0
            )
            write_row(out_csv, (grid.t, Ebase, E))
        end

        if out_every_1d > 0 && mod(step, out_every_1d) == 0
            append_data(out_h5, grid, step)
        end

        if grid.t >= next_checkpoint_t
            write_checkpoint(checkpoint_dir, grid, step, params)
            next_checkpoint_t += checkpoint_every_dt
        end

        # check if the simulation is finished
        if (max_step > 0 && step >= max_step)
            println("Max step reached, break")
            write_checkpoint(checkpoint_dir, grid, step, params)
            break
        end

        if stop_time > 0.0 && grid.t >= stop_time
            println("Stop time reached, break")
            write_checkpoint(checkpoint_dir, grid, step, params)
            break
        end

        # if slurm job time smaller than 30 mins
        if time() > slurm_job_end_time - 1800
            println("WARNING: Slurm job time is less than 30 mins, break")
            write_checkpoint(checkpoint_dir, grid, step, params)
            break
        end

        # if termination signal received
        if p.termination[]
            println("Termination signal received, break")
            write_checkpoint(checkpoint_dir, grid, step, params)
            break
        end

        apply_reflective_boundary_condition!(grid)
        step!(grid, wave_rhs!, p; mongwane=mongwane, apply_trans_zone=apply_trans_zone)

        # nan check
        if nan_check(grid)
            println("Nan detected at t = $(grid.t), step = $(step)!")
            break
        end

        step += 1
    end

    elapsed_seconds = time() - runtime_start
    elapsed = durationstring(elapsed_seconds)
    t_per_minutes = grid.t / elapsed_seconds * 60
    println(
        "Simulation finished in $elapsed, speed = $t_per_minutes, t = $(grid.t), step = $step",
    )

    ########
    # Done #
    ########
    println(
        "--------------------------------------------------------------------------------"
    )
    println("Successfully Done.")

    close(out_csv)
    close(out_h5)

    return nothing
end
