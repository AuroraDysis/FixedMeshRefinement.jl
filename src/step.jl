export step!

"""
    step!(
        grid::Grid,
        f::Function,
        p;
        mongwane::Bool=false,
        apply_trans_zone::Bool=false,
    )

Advance the solution on the entire `grid` by one time step of the coarsest level.
This function orchestrates the time stepping across all refinement levels, handling
subcycling, prolongation, restriction, and application of transition zones.

# Arguments
- `grid::Grid`: The FMR grid structure.
- `f::Function`: The function that computes the right-hand side of the ODEs.
  It should have the signature `f(level, k, u, p, t)`.
- `p`: Parameters to be passed to the RHS function `f`.
- `mongwane::Bool`: If `true`, enables Mongwane's subcycling method. Defaults to `false`.
- `apply_trans_zone::Bool`: If `true`, applies transition zones to smooth inter-grid
  boundaries. Defaults to `false`.
"""
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
                prolongate_mongwane!(grid, l, false)
            else
                prolongate!(grid, l, false)
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
        substeps = ones(Int, max_level)
        for s in 2:(2^(max_level - 1))  # from second to final substep (of the finest level)
            for l in 2:max_level  # march all levels except the coarest (from coarse to fine)
                if l == max_level ||
                    (2 * substeps[l] == substeps[l + 1] && substeps[l] < 2^(l - 1))
                    substeps[l] += 1
                    interp_in_time = iseven(substeps[l])

                    if l < max_level
                        restrict_injection!(grid, l; apply_trans_zone=apply_trans_zone)  # from l+1 to l
                    end

                    # from l-1 to l
                    if mongwane
                        prolongate_mongwane!(grid, l, interp_in_time)
                    else
                        prolongate!(grid, l, interp_in_time)
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
        restrict_injection!(grid, l; apply_trans_zone=apply_trans_zone)
    end

    #------------------#
    # update grid t #
    #------------------#
    grid.t = grid.levels[1].t

    return nothing
end
