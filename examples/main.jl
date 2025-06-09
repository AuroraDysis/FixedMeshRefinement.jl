using Infino
using Printf
using TOML

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
    diss = get(params, "diss", 0.0)
    subcycling = get(params, "subcycling", true)
    mongwane = get(params, "mongwane", false)
    ntrans = get(params, "ntrans", 3)
    ord_s = get(params, "ord_s", 3)
    apply_trans_zone = get(params, "apply_trans_zone", false)
    initial_data = get(params, "initial_data", "Gaussian")
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
    grid = Infino.Basic.Grid(
        num_interior_points,
        domain_boxes,
        num_ghost_points,
        num_buffer_points;
        ntrans=ntrans,
        ord_s=ord_s,
        cfl=cfl,
        diss=diss,
        subcycling=subcycling,
    )
    gfs = Infino.Basic.GridFunction(2, grid)

    ###############
    # Intial Data #
    ###############
    println("Setting up initial conditions...")
    println("  initial data type: $initial_data")
    if initial_data == "Gaussian"
        Infino.InitialData.Gaussian!(gfs)
    elseif initial_data == "sinusoidal"
        Infino.InitialData.sinusoidal!(gfs)
    else
        println("Initial data type '$initial_data' unsupported yet")
        exit()
    end

    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    if !Mongwane
        Infino.InitialData.MarchBackwards!(gfs)
        Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    end

    @printf(
        "Simulation time: %.4f, iteration %d. E = %.4f\n",
        grid.t,
        0,
        Infino.Physical.Energy(gfs)
    )
    Infino.WriteIO.dump(out_dir, gfs, 0)

    ##########
    # Evolve #
    ##########
    println("Start evolution...")

    for i in 1:max_step
        Infino.ODESolver.Evolve!(
            Infino.Physical.WaveRHS!,
            gfs;
            mongwane=mongwane,
            apply_trans_zone=apply_trans_zone,
        )

        @printf(
            "Simulation time: %.4f, iteration %d. E = %.4f\n",
            grid.t,
            i,
            Infino.Physical.Energy(gfs)
        )
        if (mod(i, out_every) == 0)
            Infino.WriteIO.dump(out_dir, gfs, i)
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
pars_path = ARGS[1]
params = TOML.parsefile(pars_path)

# create output directory
out_dir = joinpath(
    dirname(pars_path),
    haskey(params, "out_dir") ? params["out_dir"] : splitext(basename(pars_path))[1],
)
if isdir(out_dir)
    println("Removing old directory '$out_dir'...")
    rm(out_dir; recursive=true)
end
println("Creating new directory '$out_dir'...")
mkdir(out_dir)

# copy parfile into out_dir
cp(pars_path, out_dir * "/" * basename(pars_path))

# config
redirect_std =
    haskey(params["configs"], "redirect_std") ? params["configs"]["redirect_std"] : true

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
