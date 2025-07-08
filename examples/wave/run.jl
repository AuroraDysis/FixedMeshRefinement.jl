include("main.jl")

using JLD2
using TOML
using ThreadPinning

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

function main()
    # check arguments
    if length(ARGS) < 1
        println("Usage: julia run.jl parfile.toml")
        println("Or:    julia run.jl checkpoint.jld2")
        exit(1)
    end

    # check if the file exists
    if !isfile(ARGS[1])
        println("Error: file '$(ARGS[1])' not found.")
        exit(1)
    end

    # Pin threads
    if ThreadPinning.SLURM.hasfullnode()
        pinthreads(:cores)
    else
        pinthreads(:affinitymask)
    end

    input_file = ARGS[1]

    if endswith(input_file, ".toml")
        params_path = input_file
        params = TOML.parsefile(params_path)

        # create output directory
        out_dir = params["out_dir"]

        if !isabspath(out_dir)
            out_dir = joinpath(
                dirname(params_path),
                get(params, "out_dir", splitext(basename(params_path))[1]),
            )
        end

        # remove out_dir if it exists
        if isdir(out_dir)
            rm(out_dir, recursive=true)
        end

        if !isdir(out_dir)
            println("Creating new directory '$out_dir'...")
            mkdir(out_dir)
        end

        # if params_path is not same as out_dir, copy parfile into out_dir
        if dirname(params_path) != out_dir
            cp(params_path, joinpath(out_dir, "config.toml"))
        end

        # config
        redirect_std = get(params, "redirect_std", true)

        if redirect_std
            # redirect output and error
            redirect_to_files(out_dir * "/stdout.txt", out_dir * "/stderr.txt") do
                pde_main(params, out_dir)
            end
        else
            pde_main(params, out_dir)
        end
    elseif endswith(input_file, ".jld2")
        checkpoint_file = input_file
        println("Restarting from checkpoint: $checkpoint_file")

        grid, step, params = jldopen(checkpoint_file, "r") do file
            file["grid"], file["step"], file["params"]
        end

        # the checkpoint file is located in out_dir/checkpoints/*.jld2
        out_dir = splitdir(dirname(checkpoint_file))[1]

        redirect_std = get(params, "redirect_std", true)

        if redirect_std
            # redirect output and error
            redirect_to_files(
                out_dir * "/stdout.txt", out_dir * "/stderr.txt"; append=true
            ) do
                pde_main(params, out_dir; grid=grid, start_step=step)
            end
        else
            pde_main(params, out_dir; grid=grid, start_step=step)
        end
    else
        println(
            "Error: unsupported file extension. Use .toml for new runs or .jld2 for restarts.",
        )
        exit(1)
    end
end

main()
