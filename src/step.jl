export step!

function step!(
    grid::Grid{NumState,NumDiagnostic},
    f::Function,
    p;
    mongwane::Bool=false,
    apply_trans_zone::Bool=false,
) where {NumState,NumDiagnostic}
    max_level = length(grid.levels)

    #-------------------------------------------------#
    # march the first substep for all levels          #
    #-------------------------------------------------#
    for l in 1:max_level  # notice that we march coarse level first
        if l > 1
            if mongwane
                prolongation_mongwane!(grid, l, false)
            else
                prolongation!(grid, l, false)
            end
            if apply_trans_zone
                apply_transition_zone!(grid, l, false)
            end
        end
        rk4!(grid.levels[l], f, p; mongwane=mongwane)
    end

    #-------------------------------------------------#
    # march the other substeps to the same time slice #
    #-------------------------------------------------#
    if grid.subcycling
        levels = grid.levels
        dt_min = levels[max_level].dt
        substeps = ones(Int, max_level)
        for s in 2:(2^(max_level - 1))  # from second to final substep (of the finest level)
            for l in 2:max_level  # march all levels except the coarest (from coarse to fine)
                if l == max_level || (
                    isapprox(levels[l].t, levels[l + 1].t; rtol=1e-12) &&
                    abs(levels[l].t - levels[1].t) > dt_min
                )
                    substeps[l] += 1
                    interp_in_time = mod(substeps[l], 2) == 0

                    if l < max_level
                        restriction!(grid, l; apply_trans_zone=apply_trans_zone)  # from l+1 to l
                    end

                    # from l-1 to l
                    if mongwane
                        prolongation_mongwane!(grid, l, interp_in_time)
                    else
                        prolongation!(grid, l, interp_in_time)
                    end

                    if apply_trans_zone
                        apply_transition_zone!(grid, l, interp_in_time)
                    end

                    rk4!(grid.levels[l], f, p; mongwane=mongwane)
                end
            end
        end
    end

    #------------------------#
    # restriction all levels #
    #------------------------#
    for l in (max_level - 1):-1:1  # notice that we restrict fine level first
        restriction!(grid, l; apply_trans_zone=apply_trans_zone)
    end

    #------------------#
    # update grid t #
    #------------------#
    grid.t = grid.levels[1].t

    return nothing
end
