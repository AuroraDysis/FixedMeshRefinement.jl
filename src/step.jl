function step!(
    grid::Grid{NumState,NumDiagnostic}; mongwane=false, apply_trans_zone=false
) where {NumState,NumDiagnostic}
    max_level = length(grid.levels)

    #-------------------------------------------------#
    # march the first substep for all levels          #
    #-------------------------------------------------#
    for l in 1:max_level  # notice that we march coarse level first
        if l > 1
            if mongwane
                Sync.prolongation_mongwane(grid, l, false)
            else
                Sync.prolongation(grid, l, false)
            end
            if apply_trans_zone
                Sync.apply_transition_zone(grid, l, false)
            end
        end
        mongwane ? rk4_Mongwane!(f, grid.levels[l]) : rk4!(f, grid.levels[l])
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
                    isapprox(levels[l].time, levels[l + 1].time; rtol=1e-12) &&
                    abs(levels[l].time - levels[1].time) > dt_min
                )
                    substeps[l] += 1
                    if l < max_level
                        Sync.restriction(grid, l; apply_trans_zone=apply_trans_zone)  # from l+1 to l
                    end
                    # from l-1 to l
                    if mongwane
                        Sync.prolongation_mongwane(grid, l, mod(substeps[l], 2) == 0)
                    else
                        Sync.prolongation(grid, l, mod(substeps[l], 2) == 0)
                    end
                    if apply_trans_zone
                        Sync.apply_transition_zone(grid, l, mod(substeps[l], 2) == 0)
                    end
                    mongwane ? rk4_Mongwane!(f, grid.levels[l]) : rk4!(f, grid.levels[l])
                end
            end
        end
    end

    #------------------------#
    # restriction all levels #
    #------------------------#
    for l in (max_level - 1):-1:1  # notice that we restrict fine level first
        Sync.restriction(grid, l; apply_trans_zone=apply_trans_zone)
    end

    #------------------#
    # update grid time #
    #------------------#
    grid.time = grid.levels[1].time

    return nothing
end
