function write_output(dir_path, grid, it)
    if isdir(dir_path)
        fname = dir_path * "/u.it" * lpad(string(it), 6, '0') * ".tsv"
        open(fname, "w") do file
            println(file, "# 1:iteation 2:time 3:level 4:i 5:x 6:psi 7:Pi")

            (; num_levels, levels) = grid
            for l in 1:num_levels
                level = levels[l]
                u = level.state[end]
                x = level.x
                for i in eachindex(x)
                    println(
                        file,
                        it,
                        " ",
                        level.t,
                        " ",
                        l,
                        " ",
                        i,
                        " ",
                        x[i],
                        " ",
                        u[i, 1],
                        " ",
                        u[i, 2],
                    )
                end
            end
        end
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end
