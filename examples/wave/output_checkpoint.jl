using JLD2

function write_checkpoint(dir_path, grid, step, params)
    if isdir(dir_path)
        checkpoint_file = joinpath(dir_path, "checkpoint.step$(lpad(step, 8, '0')).t$(lpad(@sprintf("%.3f", grid.t), 9, '0')).jld2")
        jldsave(checkpoint_file; grid=grid, step=step, params=params)
    else
        error("Error: directory '$dir_path' does not exist!")
    end
end
