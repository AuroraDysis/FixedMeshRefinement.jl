using FixedMeshRefinement
using Printf
using TOML

const NumState = 2
const NumDiagnostic = 1

include("initial_data.jl")

function get(params, key, default)
    return haskey(params, key) ? params[key] : default
end

function main(params, out_dir)
    ########################
    # Read Parameter Files #
    ########################
    num_interior_points = params["num_interior_points"]
    num_ghost_points = params["num_ghost_points"]
    num_buffer_points = params["num_buffer_points"]
    stop_time = get(params, "stop_time", -1.0)
    max_step = get(params, "max_step", -1)
    out_every = params["out_every"]
    domain_boxes = params["domain_boxes"]
    cfl = get(params, "cfl", 0.25)
    dissipation = get(params, "dissipation", 0.0)
    subcycling = get(params, "subcycling", true)
    mongwane = get(params, "mongwane", false)
    num_transition_points = get(params, "num_transition_points", 3)
    spatial_interpolation_order = get(params, "spatial_interpolation_order", 5)
    apply_trans_zone = get(params, "apply_trans_zone", false)
    initial_data = get(params, "initial_data", "gaussian")
    println("Parameters:")
    println("  cfl        = ", cfl)
    println("  mongwane   = ", mongwane)
    println("  trans_zone = ", apply_trans_zone)
    println("  max_step   = ", max_step)
    println("  stop_time  = ", stop_time)
    println("  out_every  = ", out_every)
    println("  out_dir    = ", out_dir)

    ########################
    # build grid structure #
    ########################
    grid = Grid{NumState,NumDiagnostic}(
        num_interior_points,
        domain_boxes,
        num_ghost_points,
        num_buffer_points;
        num_transition_points=num_transition_points,
        spatial_interpolation_order=spatial_interpolation_order,
        cfl=cfl,
        dissipation=dissipation,
        subcycling=subcycling,
    )

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
        println("Initial data type '$initial_data' unsupported yet")
        exit()
    end

    apply_periodic_boundary_condition!(grid)
    if !mongwane
        march_backwards!(grid)
        apply_periodic_boundary_condition!(grid)
    end

    @printf("Simulation time: %.4f, iteration %d. E = %.4f\n", grid.t, 0, wave_energy(grid))
    write_output(out_dir, grid, 0)

    ##########
    # Evolve #
    ##########
    println("Start evolution...")

    for i in 1:max_step
        step!(grid, wave_rhs!; mongwane=mongwane, apply_trans_zone=apply_trans_zone)

        @printf(
            "Simulation time: %.4f, iteration %d. E = %.4f\n",
            grid.t,
            i,
            wave_energy(grid)
        )
        if (mod(i, out_every) == 0)
            write_output(out_dir, grid, i)
        end

        if stop_time > 0.0 && grid.t >= stop_time
            break
        end
    end

    ########
    # Done #
    ########
    println(
        "--------------------------------------------------------------------------------"
    )
    println("Successfully Done.")

    return nothing
end

function redirect_to_files(dofunc, outfile, errfile)
    open(outfile, "w") do out
        open(errfile, "w") do err
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
    println("Usage: julia Subcycling.jl parfile.toml")
    exit(1)
end

params_path = ARGS[1]
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
    redirect_to_files("./stdout.txt", "./stderr.txt") do
        main(params, out_dir)
    end
    mv("./stdout.txt", out_dir * "/stdout.txt")
    mv("./stderr.txt", out_dir * "/stderr.txt")
else
    main(params, out_dir)
end
#===============================================================================
End Execution
===============================================================================#
