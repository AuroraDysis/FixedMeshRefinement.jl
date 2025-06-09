using Infino
using Printf
using TOML

function main(params, out_dir)
    ########################
    # Read Parameter Files #
    ########################
    num_interior_points = params["parameters"]["num_interior_points"]
    num_ghost_points = params["parameters"]["num_ghost_points"]
    num_buffer_points = params["parameters"]["num_buffer_points"]
    itlast = params["parameters"]["itlast"]
    out_every = params["parameters"]["out_every"]
    domain_boxes = params["parameters"]["domain_boxes"]
    cfl = haskey(params["parameters"], "cfl") ? params["parameters"]["cfl"] : 0.25
    diss = haskey(params["parameters"], "diss") ? params["parameters"]["diss"] : 0.0
    subcycling =
        haskey(params["parameters"], "subcycling") ? params["parameters"]["subcycling"] : true
    Mongwane =
        haskey(params["parameters"], "Mongwane") ? params["parameters"]["Mongwane"] : false
    ntrans = haskey(params["parameters"], "ntrans") ? params["parameters"]["ntrans"] : 3
    ord_s = haskey(params["parameters"], "ord_s") ? params["parameters"]["ord_s"] : 3
    apply_trans_zone =
        haskey(params["parameters"], "apply_trans_zone") ?
        params["parameters"]["apply_trans_zone"] : false
    initial_data =
        haskey(params["parameters"], "initial_data") ? params["parameters"]["initial_data"] :
        "Gaussian"
    println("Parameters:")
    println("  cfl        = ", cfl)
    println("  Mongwane   = ", Mongwane)
    println("  trans_zone = ", apply_trans_zone)
    println("  itlast     = ", itlast)
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
        ntrans = ntrans,
        ord_s = ord_s,
        cfl = cfl,
        diss = diss,
        subcycling = subcycling,
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
        gfs.grid.time,
        0,
        Infino.Physical.Energy(gfs)
    )
    Infino.WriteIO.dump(out_dir, gfs, 0)

    ##########
    # Evolve #
    ##########
    println("Start evolution...")

    for i = 1:itlast
        Infino.ODESolver.Evolve!(
            Infino.Physical.WaveRHS!,
            gfs;
            Mongwane = Mongwane,
            apply_trans_zone = apply_trans_zone,
        )

        @printf(
            "Simulation time: %.4f, iteration %d. E = %.4f\n",
            gfs.grid.time,
            i,
            Infino.Physical.Energy(gfs)
        )
        if (mod(i, out_every) == 0)
            Infino.WriteIO.dump(out_dir, gfs, i)
        end
    end

    ########
    # Done #
    ########
    println(
        "--------------------------------------------------------------------------------",
    )
    println("Successfully Done.")
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
    haskey(params["parameters"], "out_dir") ? params["parameters"]["out_dir"] :
    splitext(basename(pars_path))[1],
)
if isdir(out_dir)
    println("Removing old directory '$out_dir'...")
    rm(out_dir, recursive = true)
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
