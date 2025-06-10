using HDF5

function write_output(dir_path, grid, it)
    if isdir(dir_path)
        (; num_levels, levels) = grid

        h5open(joinpath(dir_path, "data.it$(lpad(it, 6, '0')).h5"), "w") do file
            for l in 1:num_levels
                level = levels[l]
                u = level.state[end]
                x = level.x

                g = create_group(file, "level$(lpad(l, 2, '0'))")

                write(g, "t", level.t)
                write(g, "x", collect(x))
                write(g, "psi", @view(u[:, 1]))
                write(g, "Pi", @view(u[:, 2]))
            end
        end
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end
