using HDF5
using JLD2

function write_output(dir_path, grid, it)
    if isdir(dir_path)
        (; num_levels, levels) = grid

        h5open(joinpath(dir_path, "data.it$(lpad(it, 6, '0')).h5"), "w") do file
            for l in 1:num_levels
                level = levels[l]
                x = get_x(level)
                u = get_state(level)
                num_interior_indices = get_interior_indices(level)

                g = create_group(file, "level$(lpad(l, 2, '0'))")

                write(g, "t", level.t)
                write(g, "x", collect(x[num_interior_indices]))
                write(g, "psi", @view(u[num_interior_indices, 1]))
                write(g, "Pi", @view(u[num_interior_indices, 2]))
            end
        end
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end

function write_checkpoint(dir_path, grid, step, params)
    if isdir(dir_path)
        checkpoint_file = joinpath(dir_path, "checkpoint.it$(lpad(step, 6, '0')).jld2")
        jldsave(checkpoint_file; grid=grid, step=step, params=params)
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end
