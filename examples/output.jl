module WriteIO

function write_output(dir_path, grid, it)
    if isdir(dir_path)
        fname = dir_path * "/u.it" * lpad(string(it), 6, '0') * ".tsv"
        open(fname, "w") do file
            println(file, "# 1:iteation 2:time 3:level 4:i 5:x 6:psi 7:Pi")
            for l in 1:length(grid.levels)
                level = grid.levels[l]
                levfs = grid.levels[l]
                for i in 1:length(levfs.x)
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
                        levfs.x[i],
                        " ",
                        levfs.u[1][i],
                        " ",
                        levfs.u[2][i],
                    )
                end
            end
        end
    else
        println("Error: directory '$dir_path' does not exist!")
        exit()
    end
end

end
