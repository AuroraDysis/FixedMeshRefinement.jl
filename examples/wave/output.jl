using HDF5

function write_output(dir_path, grid, it)
    if isdir(dir_path)
        (; num_levels, levels) = grid

        h5open(joinpath(dir_path, "data.it$(lpad(it, 6, '0')).h5"), "w") do file
            for l in 1:num_levels
                level = levels[l]
                u = level_state(level)
                x = level_x(level)
                num_interior_indices = level_interior_indices(level)

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
