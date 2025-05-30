export amr_main

"""
shift_field!(grid, field_index, time_levels)
Move time history one slot deeper so that the newest value can be written
into the 0-th slot.  *grid.grid_funcs* is arranged as
[field_index + i] for history i.
"""
function shift_field!(grid::AMRGrid, field_index::Int, time_levels::Int)
    # Work in Julia indices (offset +1)
    start_idx = field_index + 1
    for i in (start_idx+time_levels-1):-1:(start_idx+1)
        grid.grid_funcs[i] .= grid.grid_funcs[i-1]
    end
    return nothing
end

"""
shift_fields_one_time_level!(fields, grid)
Shift *all* fields on a single grid.
"""
function shift_fields_one_time_level!(fields::AMRField, grid::AMRGrid)
    fld = fields
    while fld !== nothing
        shift_field!(grid, fld.index, fld.time_levels)
        fld = fld.next
    end
    return nothing
end

"""
shift_grids_one_time_level!(gh)
Apply `shift_fields_one_time_level!` to every grid in the hierarchy.
"""
function shift_grids_one_time_level!(gh::AMRGridHierarchy)
    grid = gh.grids
    while grid !== nothing
        shift_fields_one_time_level!(gh.fields, grid)
        grid = grid.child
    end
    return nothing
end

# ----------------------------------------------------------------------------------
# Boundary handling for hyperbolic variables ---------------------------------------
# ----------------------------------------------------------------------------------
function quadratic_interp(parent_vals::Vector{Float64}, idx::Int, tC::Int)
    c0 = parent_vals[idx+2]                 # index+1 in C is idx+1; shift by +1 for Julia
    c1 = (parent_vals[idx+1] - parent_vals[idx+3]) / (2 * REFINE)
    c2 = (parent_vals[idx+1] - 2 * parent_vals[idx+2] + parent_vals[idx+3]) /
         (2 * REFINE^2)
    return (tstep) -> c0 + c1 * (tstep) + c2 * (tstep - 1)^2  # convenience closure
end

function set_interior_hyperbolic_boundary!(field::AMRField, parent::AMRGrid, grid::AMRGrid)
    idx = field.index + 1
    # left boundary ----------------------------------------------------------
    if grid.is_physical_boundary[1]
        perim = grid.perim_coords[1] + 1         # convert to 1-based
        tC = mod(grid.tC, REFINE)
        interp = quadratic_interp(parent.grid_funcs[idx], perim - 1, tC)
        grid.grid_funcs[idx][1] = interp(tC + 1)
        grid.grid_funcs[idx+1][1] = interp(tC)
    end
    # right boundary ---------------------------------------------------------
    if grid.is_physical_boundary[2]
        perim = grid.perim_coords[2] + 1
        tC = mod(grid.tC, REFINE)
        interp = quadratic_interp(parent.grid_funcs[idx], perim - 1, tC)
        last = grid.Nx
        grid.grid_funcs[idx][last] = interp(tC + 1)
        grid.grid_funcs[idx+1][last] = interp(tC)
    end
    return nothing
end

function apply_interior_hyperbolic_boundaries!(fields::AMRField, parent::AMRGrid, grid::AMRGrid)
    fld = fields
    while fld !== nothing
        if fld.pde_type == HYPERBOLIC
            set_interior_hyperbolic_boundary!(fld, parent, grid)
        end
        fld = fld.next
    end
    return nothing
end

# ----------------------------------------------------------------------------------
# ODE extrapolation helpers --------------------------------------------------------
# ----------------------------------------------------------------------------------
function linear_extrap_level0!(field::AMRField, grid::AMRGrid)
    idx = field.index + 1
    extrap_idx = idx + field.time_levels
    for j in 1:grid.Nx
        p0 = grid.grid_funcs[extrap_idx][j]
        p1 = p0 - grid.grid_funcs[extrap_idx+1][j]
        grid.grid_funcs[idx][j] = p0 + p1
    end
    return nothing
end

function linear_extrap_leveln!(field::AMRField, grid::AMRGrid)
    idx = field.index + 1
    extrap_idx = idx + field.time_levels
    step = (grid.tC == grid.parent.tC * REFINE) ? REFINE : mod(grid.tC, REFINE)
    for j in 1:grid.Nx
        p0 = grid.grid_funcs[extrap_idx][j]
        p1 = (p0 - grid.grid_funcs[extrap_idx+1][j]) / REFINE
        grid.grid_funcs[idx][j] = p0 + p1 * step
    end
    return nothing
end

function amr_extrapolate_ode_fields!(fields::AMRField, grid::AMRGrid)
    fld = fields
    while fld !== nothing
        if fld.pde_type == ODE
            if grid.level == 0
                linear_extrap_level0!(fld, grid)
            else
                linear_extrap_leveln!(fld, grid)
            end
        end
        fld = fld.next
    end
    return nothing
end

# ----------------------------------------------------------------------------------
# Higher-level ODE solver orchestration -------------------------------------------
# ----------------------------------------------------------------------------------
function set_finer_grid_ode_ic!(fields::AMRField, parent::AMRGrid, child::AMRGrid)
    fld = fields
    child_lower = child.perim_coords[1] + 1
    while fld !== nothing
        if fld.pde_type == ODE
            idx = fld.index + 1
            child.grid_funcs[idx][1] = parent.grid_funcs[idx][child_lower]
        end
        fld = fld.next
    end
end

function set_coarser_grid_ode_ic!(fields::AMRField, parent::AMRGrid, child::AMRGrid)
    fld = fields
    child_upper = child.perim_coords[2] + 1
    while fld !== nothing
        if fld.pde_type == ODE
            idx = fld.index + 1
            parent.grid_funcs[idx][child_upper] = child.grid_funcs[idx][child.Nx]
        end
        fld = fld.next
    end
end

# Recursive ODE solving along hierarchy ------------------------------------------
function solve_ode_fields!(fields::AMRField, grid::AMRGrid, solve_ode_fn::Function)
    # If finest grid => just solve and return
    if grid.child === nothing
        solve_ode_fn(grid)
        return
    end

    excised = grid.excised_jC
    child_lower = grid.child.perim_coords[1]
    child_upper = grid.child.perim_coords[2]

    if child_upper < excised
        solve_ode_fn(grid)
        return
    end

    # Solve portion to the left of child
    if child_lower > excised
        original_bbox = grid.bbox
        original_Nx = grid.Nx

        grid.bbox = (grid.bbox[1], grid.bbox[1] + child_lower * grid.dx)
        grid.Nx = child_lower + 1
        solve_ode_fn(grid)

        grid.bbox = original_bbox
        grid.Nx = original_Nx
    end

    # Solve child (ensure IC)
    if child_lower > 0
        set_finer_grid_ode_ic!(fields, grid, grid.child)
        solve_ode_fields!(fields, grid.child, solve_ode_fn)
    end

    # Solve RHS of child
    if child_upper < grid.Nx - 1
        grid.excised_jC = child_upper
        grid.excision_on = false
        set_coarser_grid_ode_ic!(fields, grid, grid.child)
        solve_ode_fn(grid)
        grid.excised_jC = excised
        grid.excision_on = true
    end
end

# -------------------------------------------------------------------------------
# Time-integration driver (recursive across levels) -----------------------------
# -------------------------------------------------------------------------------
function evolve_grid!(fields::AMRField, grid::AMRGrid, num_t_steps::Int,
    evolve_hyp_fn::Function, solve_ode_fn::Function)
    for _ in 1:num_t_steps
        # Interior boundary conditions for hyperbolics
        if (grid.parent !== nothing) && (grid.level > 1)
            apply_interior_hyperbolic_boundaries!(fields, grid.parent, grid)
        end

        # Advance time counters
        grid.tC += 1
        grid.time += grid.dt

        shift_fields_one_time_level!(fields, grid)
        amr_extrapolate_ode_fields!(fields, grid)

        evolve_hyp_fn(grid)                      # user-supplied hyperbolic solver

        if grid.child !== nothing
            evolve_grid!(fields, grid.child, REFINE, evolve_hyp_fn, solve_ode_fn)
        end

        if grid.level > 0
            solve_ode_fields!(fields, grid, solve_ode_fn)
        end
    end

    # Inject synchronised grids
    if grid.child !== nothing
        inject_overlaping_fields(fields, grid.child, grid)
    end

    # TODO: set_grid_ode_extrap_levels! (omitted for brevity)
    return nothing
end

# -------------------------------------------------------------------------------
# Public interface --------------------------------------------------------------
# -------------------------------------------------------------------------------
function amr_main(gh::AMRGridHierarchy,
    free_initial_data::Function,
    evolve_hyperbolic_pde::Function,
    solve_ode::Function,
    compute_diagnostics::Function,
    save_to_file::Function)

    # (1) Add initial nested grids (fixed layout for now)
    add_initial_grids(gh)

    # (2) Placeholder: user initial data on all grids (simplified)
    grid = gh.grids
    while grid !== nothing
        free_initial_data(grid)
        grid = grid.child
    end

    # (3) Save initial state
    save_to_file(gh.grids)

    # (4) Main evolution loop --------------------------------------------------
    for step in 1:gh.Nt
        evolve_grid!(gh.fields, gh.grids, 1, evolve_hyperbolic_pde, solve_ode)

        if step % gh.t_step_save == 0
            compute_diagnostics(gh.grids)
            save_to_file(gh.grids)
        end
    end
end
