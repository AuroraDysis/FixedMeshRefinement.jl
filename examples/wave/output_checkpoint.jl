using JLD2

function write_checkpoint(dir_path, grid, step, params)
    if isdir(dir_path)
        checkpoint_file = joinpath(dir_path, "checkpoint.it$(lpad(step, 6, '0')).jld2")
        jldsave(checkpoint_file; grid=grid, step=step, params=params)
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end
