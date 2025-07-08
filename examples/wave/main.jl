include("../../src/FixedMeshRefinement.jl")
using .FixedMeshRefinement

using Printf
using TOML
using JLD2
using FastBroadcast

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

function main(params, out_dir; grid=nothing, start_step=1)
    ########################
    # Read Parameter Files #
    ########################
    stop_time = get(params, "stop_time", -1.0)
    max_step = get(params, "max_step", -1)
    out_every_0d = get(params, "out_every_0d", 10)
    out_every_1d = get(params, "out_every_1d", 100)
    cfl = get(params, "cfl", 0.25)
    dissipation = get(params, "dissipation", 0.0)
    subcycling = get(params, "subcycling", true)
    mongwane = get(params, "mongwane", false)
    apply_trans_zone = get(params, "apply_trans_zone", true)
    checkpoint_every = get(params, "checkpoint_every", -1)

    println("Parameters:")
    println("  cfl              = ", cfl)
    println("  mongwane         = ", mongwane)
    println("  trans_zone       = ", apply_trans_zone)
    println("  max_step         = ", max_step)
    println("  stop_time        = ", stop_time)
    println("  out_every_0d     = ", out_every_0d)
    println("  out_every_1d     = ", out_every_1d)
    println("  checkpoint_every = ", checkpoint_every)
    println("  out_dir          = ", out_dir)

    ########################
    # build grid structure #
    ########################
    p = (; dissipation)
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

    step = start_step
    while (max_step > 0 && step <= max_step) || (stop_time > 0.0 && grid.t < stop_time)
        apply_reflective_boundary_condition!(grid)
        step!(grid, wave_rhs!, p; mongwane=mongwane, apply_trans_zone=apply_trans_zone)

        Ebase, E = wave_energy(grid)
        @printf(
            "t = %.4f, iteration %d. dEbase = %.5g, dE = %.5g\n",
            grid.t,
            step,
            Ebase - E0,
            E - E0
        )

        if out_every_0d > 0 && mod(step, out_every_0d) == 0
            write_row(out_csv, (grid.t, Ebase, E))
        end

        if out_every_1d > 0 && mod(step, out_every_1d) == 0
            append_data(out_h5, grid, step)
        end

        if checkpoint_every > 0 && mod(step, checkpoint_every) == 0
            write_checkpoint(out_dir, grid, step, params)
        end

        # nan check
        if nan_check(grid)
            println("Nan detected at t = $(grid.t), step = $(step)!")
            exit()
        end

        step += 1
    end

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

function redirect_to_files(dofunc, outfile, errfile; append=false)
    mode = append ? "a" : "w"
    open(outfile, mode) do out
        open(errfile, mode) do err
            redirect_stdout(out) do
                redirect_stderr(err) do
                    dofunc()
                end
            end
        end
    end
end

#===============================================================================
Start Execution
===============================================================================#
if length(ARGS) < 1
    println("Usage: julia main.jl parfile.toml")
    println("Or:    julia main.jl checkpoint.jld2")
    exit(1)
end

input_file = ARGS[1]
if endswith(input_file, ".toml")
    params_path = input_file
    params = TOML.parsefile(params_path)

    # create output directory
    out_dir = joinpath(
        dirname(params_path), get(params, "out_dir", splitext(basename(params_path))[1])
    )
    if isdir(out_dir)
        println("Removing old directory '$out_dir'...")
        rm(out_dir; recursive=true)
    end
    println("Creating new directory '$out_dir'...")
    mkdir(out_dir)

    # copy parfile into out_dir
    cp(params_path, out_dir * "/" * basename(params_path))

    # config
    redirect_std = get(params, "redirect_std", true)

    if redirect_std
        # redirect output and error
        redirect_to_files(out_dir * "/stdout.txt", out_dir * "/stderr.txt") do
            main(params, out_dir)
        end
    else
        main(params, out_dir)
    end
elseif endswith(input_file, ".jld2")
    checkpoint_file = input_file
    println("Restarting from checkpoint: $checkpoint_file")
    grid, step, params = jldopen(checkpoint_file, "r") do file
        file["grid"], file["step"], file["params"]
    end

    out_dir = dirname(checkpoint_file)
    redirect_std = get(params, "redirect_std", true)

    if redirect_std
        # redirect output and error
        redirect_to_files(out_dir * "/stdout.txt", out_dir * "/stderr.txt"; append=true) do
            main(params, out_dir; grid=grid, start_step=step + 1)
        end
    else
        main(params, out_dir; grid=grid, start_step=step + 1)
    end
else
    println(
        "Error: unsupported file extension. Use .toml for new runs or .jld2 for restarts."
    )
    exit(1)
end
#===============================================================================
End Execution
===============================================================================#
