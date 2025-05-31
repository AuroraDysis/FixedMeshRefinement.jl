using ArgCheck

struct AMRContext{TF<:AbstractFloat}
    # Grid layout parameters
    max_levels::Int
    refinement_ratio::Int
    regrid_cadence::Int
    buffer_coord::Int
    min_grid_size::Int

    # Simulation parameters
    fields::Vector{AMRField}

    cfl::TF # CFL number

    num_grid_functions::Int
    grid_functions_storage_indices::Vector{Int}

    # Error tolerance parameters
    trunc_err_tolerance::TF
    min_shift_distance::Int

    function AMRContext(
        fields::Vector{AMRField};
        cfl::TF=0.25,
        use_excision::Bool=false,
        max_levels::Int=8,
        refinement_ratio::Int=2,
        regrid_cadence::Int=80,
        buffer_coord::Int=40,
        min_grid_size::Int=40,
        trunc_err_tolerance::TF=1e-6,
        min_shift_distance::Int=8,
    ) where {TF<:AbstractFloat}
        num_grid_functions = [
            field.num_time_levels + field.num_extrapolation_levels for field in fields
        ]
        grid_functions_storage_indices = Vector{Int}(undef, num_grid_functions)
        current_index = 1
        for (i, field) in enumerate(fields)
            grid_functions_storage_indices[i] = current_index
            current_index += field.num_time_levels + field.num_extrapolation_levels
        end
        return new{TF}(
            max_levels,
            refinement_ratio,
            regrid_cadence,
            buffer_coord,
            min_grid_size,
            fields,
            cfl,
            num_grid_functions,
            grid_functions_storage_indices,
            trunc_err_tolerance,
            min_shift_distance,
        )
    end
end

@enum PDEType begin
    HYPERBOLIC
    ELLIPTIC
    ODE
    DIAGNOSTIC
end

mutable struct AMRField
    const name::String
    const num_time_levels::Int                 # hyperbolic time levels stored on each grid
    const num_extrapolation_levels::Int        # extra history for ODE extrapolation
    const pde_type::PDEType                    # PDE type (hyperbolic, elliptic, ODE, diagnostic)
    regrid_indices::NTuple{2,Int}         # (lower, upper) flagged coordinates for regridding

    function AMRField(
        name::String, pde_type::PDEType, num_time_levels::Int, num_extrapolation_levels::Int
    )
        return new(name, num_time_levels, num_extrapolation_levels, pde_type, [0, 0])
    end
end

mutable struct AMRLevel{TF<:AbstractFloat}
    ctx::AMRContext{TF}
    parent::Union{AMRLevel{TF},Nothing}
    child::Union{AMRLevel{TF},Nothing}
    level::Int                                                         # level of the grid
    bounding_box::NTuple{2,TF}                                         # physical domain covered by this grid
    num_grid_points::Int                                               # number of spatial grid points
    is_physical_boundary::NTuple{2,Bool}                               # whether the grid touches a physical boundary
    parent_indices::NTuple{2,Int}                                      # coordinates inside parent grid (inclusive)
    regrid_indices::NTuple{2,Int}                                      # tagged regions for regridding
    excision_index::Int                                                # excision coordinate (0 => none)
    grid_functions_storage::Vector{Vector{TF}}                         # storage: each entry is a 1-D array over space
    time::TF                                                           # current simulation time
    dt::TF                                                             # local time step
    dx::TF                                                             # spatial resolution
end

"""
    AMRLevel(ctx::AMRContext{TF}, bounding_box::NTuple{2,TF}, num_grid_points::Int) where {TF<:AbstractFloat}
Construct the base AMRLevel, with a given bounding box and number of grid points.
"""
function AMRLevel(
    ctx::AMRContext{TF}, bounding_box::NTuple{2,TF}, num_grid_points::Int
) where {TF<:AbstractFloat}
    level = 1
    parent = nothing
    child = nothing
    is_physical_boundary = (true, true)
    parent_indices = (0, 0)
    regrid_indices = (0, 0)
    excision_index = 0
    grid_functions_storage = [
        Vector{TF}(NaN, num_grid_points) for _ in 1:(ctx.num_grid_functions)
    ]
    time = zero(TF)
    dx = (bounding_box[2] - bounding_box[1]) / (num_grid_points - 1)
    dt = ctx.cfl * dx
    return new{TF}(
        ctx,
        parent,
        child,
        level,
        bounding_box,
        num_grid_points,
        is_physical_boundary,
        parent_indices,
        regrid_indices,
        excision_index,
        grid_functions_storage,
        time,
        dt,
        dx,
    )
end

"""
    AMRLevel(parent, start_idx, end_idx) -> AMRLevel
Construct a new AMRLevel refined with respect to `parent`.
The new grid covers the interval [`start_idx`, `end_idx`] in parent coordinates.
"""
function AMRLevel(
    parent::AMRLevel{TF}, start_idx::Int, end_idx::Int
) where {TF<:AbstractFloat}
    ctx = parent.ctx

    @argcheck parent.level + 1 <= ctx.max_levels "Cannot create level $(parent.level + 1), maximum AMR levels is $(ctx.max_levels)"
    @argcheck start_idx >= 1 && start_idx <= parent.num_grid_points "start_idx ($start_idx) is out of bounds [1, $(parent.num_grid_points)] for parent grid."
    @argcheck end_idx >= 1 && end_idx <= parent.num_grid_points "end_idx ($end_idx) is out of bounds [1, $(parent.num_grid_points)] for parent grid."
    @argcheck start_idx <= end_idx "start_idx ($start_idx) cannot be greater than end_idx ($end_idx)"

    # Level of the new (child) grid
    child_level_val = parent.level + 1

    # Number of spatial points in the new grid
    num_grid_points = ctx.refinement_ratio * (end_idx - start_idx) + 1

    # Spatial and temporal resolution for the new grid
    dx_val = parent.dx / ctx.refinement_ratio
    dt_val = parent.dt / ctx.refinement_ratio

    # Child grid starts at the same simulation time as parent
    time_val = parent.time

    # Physical bounding box of the new grid
    child_bbox_start = parent.bounding_box[1] + (start_idx - 1) * parent.dx
    child_bbox_end = parent.bounding_box[1] + (end_idx - 1) * parent.dx
    bounding_box_val = (child_bbox_start, child_bbox_end)

    # Parent grid bounds (inclusive indices on parent grid that this child covers)
    parent_indices_val = (start_idx, end_idx)

    # Determine if boundaries of the new grid are physical boundaries
    is_left_physical = parent.is_physical_boundary[1] && (start_idx == 1)
    is_right_physical =
        parent.is_physical_boundary[2] && (end_idx == parent.num_grid_points)
    is_physical_boundary_val = (is_left_physical, is_right_physical)

    # Excision index mapping
    if parent.excision_index >= start_idx && parent.excision_index <= end_idx
        excision_index_val = ctx.refinement_ratio * (parent.excision_index - start_idx) + 1
    else
        excision_index_val = 0
    end

    grid_functions_storage = [
        Vector{TF}(NaN, num_grid_points) for _ in 1:(ctx.num_grid_functions)
    ]

    # Create the new AMRLevel instance using the `new` keyword for constructors
    return new{TF}(
        ctx,
        parent,                             # parent
        nothing,                            # child (initialized to nothing)
        child_level_val,                    # level
        bounding_box_val,                   # bounding_box
        num_grid_points,                    # num_grid_points
        is_physical_boundary_val,           # is_physical_boundary
        parent_indices_val,                 # parent_indices
        (1, 0),                             # regrid_indices
        excision_index_val,                 # excision_index
        grid_functions_storage,             # grid_functions_storage
        time_val,                           # time
        dt_val,                             # dt
        dx_val,                             # dx
    )
end

function amr_refine_level!(parent::AMRLevel, start_idx::Int, end_idx::Int)
    parent.child = AMRLevel(parent, start_idx, end_idx)
    parent.child.parent = parent
    return nothing
end

"""
amr_insert_level!(parent, new_child_grid)
Insert `new_child_grid` directly beneath `parent`, updating the linked structure.
"""
function amr_insert_level!(parent::AMRLevel, new_child_grid::AMRLevel)
    new_child_grid.child = parent.child
    new_child_grid.parent = parent
    if !isnothing(parent.child)
        parent.child.parent = new_child_grid
    end
    parent.child = new_child_grid

    # increase the level of the new child grid
    child_of_parent = new_child_grid
    while !isnothing(child_of_parent.child)
        child_of_parent.child.level = child_of_parent.level + 1
        child_of_parent = child_of_parent.child
    end

    return nothing
end

"""
amr_restrict_level!(grid::AMRLevel)
Restrict overlapping grid-functions from `grid` into its `parent`.
This version assumes the 'index' corresponded to the first data slot of each field.
"""
function amr_restrict_to_parent!(child::AMRLevel)
    parent = child.parent

    if isnothing(parent)
        return nothing
    end

    (; ctx, parent_indices) = parent
    (; refinement_ratio, fields, grid_functions_storage_indices) = ctx

    # restrict the data from the child grid into the parent grid along shared grid points
    for i in 1:length(fields)
        storage_index = grid_functions_storage_indices[i]
        parent_data = parent.grid_functions_storage[storage_index]
        child_data = child.grid_functions_storage[storage_index]
        parent_data .= @view child_data[parent_indices[1]:refinement_ratio:parent_indices[2]]
    end

    return nothing
end

"""
amr_prolong_level!(grid::AMRLevel)
Prolong the grid-functions from `grid` into its `child`.
"""
function amr_prolong_to_child!(parent::AMRLevel)
    child = parent.child

    if isnothing(child)
        return nothing
    end

    (; ctx, parent_indices) = child
    (; refinement_ratio, num_grid_functions) = ctx

    # prolong all grid functions from parent to child
    start_idx = parent_indices[1]
    end_idx = parent_indices[2]
    num_points = end_idx - start_idx + 1
    for i in 1:num_grid_functions
        parent_data = parent.grid_functions_storage[i]
        child_data = child.grid_functions_storage[i]
        for j in 1:(num_points - 1)
            p0 = parent_data[start_idx + j - 1]
            p1 =
                (parent_data[start_idx + j] - parent_data[start_idx + j - 1]) /
                refinement_ratio
            for k in 1:refinement_ratio
                child_data[refinement_ratio * (j - 1) + k] = p0 + (k - 1) * p1
            end
        end
        child_data[refinement_ratio * (num_points - 1) + 1] = parent_data[end_idx]
    end

    return nothing
end

"""
    get_base_level(level::AMRLevel)
Navigate to the base (coarsest) level from the given `level`.
"""
function get_base_level(level::AMRLevel)
    current = level
    while !isnothing(current.parent)
        current = current.parent
    end
    return current
end

"""
    get_finest_level(level::AMRLevel)
Navigate to the finest level from the given `level`.
"""
function get_finest_level(level::AMRLevel)
    current = level
    while !isnothing(current.child)
        current = current.child
    end
    return current
end

"""
    find_level(start_level::AMRLevel, target_level_num::Int)
Find and return the level at `target_level_num`, starting the search from `start_level`.
Returns `nothing` if the level is not found.
"""
function find_level(start_level::AMRLevel, target_level_num::Int)
    @argcheck target_level_num >= 1 "Target level must be positive."

    # Navigate to base grid to ensure consistent starting point if needed,
    # or decide search direction based on current vs target.
    current = get_base_level(start_level)

    while !isnothing(current) && current.level != target_level_num
        if current.level < target_level_num
            current = current.child
        else # current.level > target_level_num, should not happen if starting from base
            return nothing # Or error, depending on desired behavior
        end
    end

    return current
end

"""
    amr_remove_level!(grid_to_remove::AMRLevel)
Remove `grid_to_remove` from the AMR hierarchy, updating links.
The children of the removed grid are not automatically handled (i.e., they are orphaned).
Consider explicit handling of children if necessary before calling this.
"""
function amr_remove_level!(grid_to_remove::AMRLevel)
    parent = grid_to_remove.parent
    child_grid = grid_to_remove.child

    if !isnothing(parent)
        parent.child = child_grid # Link parent to grandchild
    end

    if !isnothing(child_grid)
        child_grid.parent = parent # Link grandchild to grandparent
        # Update levels of subsequent child grids
        current_child = child_grid
        while !isnothing(current_child)
            if !isnothing(current_child.parent)
                current_child.level = current_child.parent.level + 1
            else
                # This case implies child_grid became the new base grid
                # This logic might need refinement if removing the base grid (level 1)
                # or if a more complex re-leveling strategy is needed.
                # For now, assume base grid is not removed or handled separately.
                current_child.level = 1 # Or adjust based on new parent.
            end
            current_child = current_child.child
        end
    end

    # Clear links from the removed grid
    grid_to_remove.parent = nothing
    grid_to_remove.child = nothing

    # The grid_to_remove is now unlinked and will be garbage collected by Julia.
    return nothing
end

# Stubs for more complex functions from amr_grid_hierarchy.c

"""
STUB: smooth_all_grid_funcs!(fields::Vector{AMRField}, grid::AMRLevel)
Apply a Kreiss-Oliger filter to the grid functions.
"""
function smooth_all_grid_funcs!(fields::Vector{AMRField}, grid::AMRLevel)
    epsilon_ko = 1.0
    Nx = grid.num_grid_points

    # Determine the actual start index of valid data (1-based)
    # If grid.excision_index is 0, it means no excision, so data starts at 1.
    # Otherwise, data starts at grid.excision_index.
    first_valid_idx = grid.excision_index == 0 ? 1 : grid.excision_index

    # Iterate over each field provided
    for (field_enum_idx, _field_obj) in enumerate(fields)
        # Get the index for the grid_functions_storage array for this field's primary data
        storage_idx = grid.ctx.grid_functions_storage_indices[field_enum_idx]
        vals = grid.grid_functions_storage[storage_idx]

        # KO filter requires at least 5 points for full stencil application.
        # The C code does not explicitly check this, implying Nx is assumed large enough.
        # We'll add checks for stencil validity.

        # Main loop for interior points
        # Loop runs if first_valid_idx + 2 <= Nx - 2  => Nx >= first_valid_idx + 4
        loop_start_interior = first_valid_idx + 2
        loop_end_interior = Nx - 2
        if loop_start_interior <= loop_end_interior
            for j in loop_start_interior:loop_end_interior
                # Ensure stencil points are within bounds [first_valid_idx, Nx]
                # For this loop, j-2 >= first_valid_idx and j+2 <= Nx is guaranteed by loop bounds
                if (j - 2 >= first_valid_idx && j + 2 <= Nx) # Redundant given loop_start/end but safe
                    stencil_val =
                        vals[j + 2] - 4 * vals[j + 1] + 6 * vals[j] - 4 * vals[j - 1] +
                        vals[j - 2]
                    vals[j] -= (epsilon_ko / 16.0) * stencil_val
                end
            end
        end

        # Boundary condition at Nx-1 (C index Nx-2)
        # Stencil: vals[Nx-4] to vals[Nx]. Target: vals[Nx-1].
        # Requires: Nx >= 5 (for 5 points).
        # And all stencil points vals[Nx-4]...vals[Nx] must be valid.
        # Smallest index used is Nx-4. Target index is Nx-1.
        if Nx >= 5 && (Nx - 1 >= first_valid_idx) && (Nx - 4 >= first_valid_idx)
            stencil_val_hi =
                vals[Nx] - 4 * vals[Nx - 1] + 6 * vals[Nx - 2] - 4 * vals[Nx - 3] +
                vals[Nx - 4]
            vals[Nx - 1] += (epsilon_ko / 16.0) * stencil_val_hi
        end

        # Boundary condition at 2 (C index 1)
        # Stencil: vals[1] to vals[5]. Target: vals[2].
        # Requires: Nx >= 5.
        # All stencil points vals[1]...vals[5] must be valid.
        # Smallest index used is 1. Target index is 2.
        if Nx >= 5 && (2 >= first_valid_idx) && (1 >= first_valid_idx) # (1 >= first_valid_idx implies first_valid_idx is 1)
            stencil_val_lo = vals[1] - 4 * vals[2] + 6 * vals[3] - 4 * vals[4] + vals[5]
            vals[2] += (epsilon_ko / 16.0) * stencil_val_lo
        end
    end
    return nothing
end

"""
STUB: tag_regridding_regions_by_richardson!(parent::AMRLevel, child::AMRLevel)
Tag regions for regridding based on Richardson extrapolation error.
Updates `regrid_indices` in `fields`.
"""
function tag_regridding_regions_by_richardson!(parent::AMRLevel, child::AMRLevel)
    (; ctx) = child
    (; fields, refinement_ratio, buffer_coord, grid_functions_storage_indices) = ctx

    for (field_idx, field) in enumerate(fields)
        if field.pde_type != HYPERBOLIC
            # For non-hyperbolic fields, we don't flag any regions.
            field.regrid_indices = (1, 0)
            continue
        end

        storage_idx = grid_functions_storage_indices[field_idx]
        parent_data = parent.grid_functions_storage[storage_idx]
        child_data = child.grid_functions_storage[storage_idx]

        # Initialize to indicate no points flagged yet.
        lower_flagged_child_coord = child.num_grid_points + 1
        upper_flagged_child_coord = 0

        # These are the indices on the parent grid that the child grid `grid` covers.
        parent_idx_start_of_child = child.parent_indices[1]
        parent_idx_end_of_child = child.parent_indices[2]

        # Determine the iteration range on the parent grid, applying buffers.
        # The C code's `grid->perim_interior` refers to the child grid's boundary type.
        # `grid.is_physical_boundary[1]` is true if child's left is physical.
        # If false, it's an interior boundary, so buffer is applied.

        iter_parent_start = parent_idx_start_of_child
        if !child.is_physical_boundary[1] # if child's left boundary is not physical (i.e., interior)
            iter_parent_start += buffer_coord
        end

        iter_parent_end = parent_idx_end_of_child
        if !child.is_physical_boundary[2] # if child's right boundary is not physical (i.e., interior)
            iter_parent_end -= buffer_coord
        end

        # C loop is `for (jC=start_jC; jC<end_jC; jC++)`, so end_jC is exclusive.
        # We iterate j_parent from iter_parent_start up to iter_parent_end - 1.
        # However, the Richardson error compares parent[j_p] with child[corresponding to j_p].
        # The loop in C is `for (int jC=start_jC; jC<end_jC; jC++)`.
        # It seems indices should cover the comparable region.
        # If iter_parent_end becomes less than iter_parent_start, loop won't run.

        for j_parent in iter_parent_start:(iter_parent_end - 1) # Adjusted for C's exclusive end
            if j_parent < 1 || j_parent > parent.num_grid_points
                # Should not happen if parent_indices are valid for parent
                continue
            end

            parent_val = parent_data[j_parent]

            # Convert parent index j_parent to child index
            # C: grid_index = refinement*(jC-lower_jC)
            # lower_jC is parent_idx_start_of_child (0-based in C, 1-based in Julia)
            # So, child_coord = refinement_ratio * (j_parent - parent_idx_start_of_child) + 1 (for 1-based child array)
            child_coord = refinement_ratio * (j_parent - parent_idx_start_of_child) + 1

            if child_coord < 1 || child_coord > grid.num_grid_points
                # This indicates an issue with index mapping or loop bounds.
                # Given the C logic, this comparison point should be valid.
                continue
            end
            child_val = child_data[child_coord]

            trunc_err = abs(parent_val - child_val)

            if trunc_err > child.ctx.trunc_err_tolerance
                if lower_flagged_child_coord > child.num_grid_points # First time flagging
                    lower_flagged_child_coord = child_coord
                    upper_flagged_child_coord = child_coord
                else
                    lower_flagged_child_coord = min(lower_flagged_child_coord, child_coord)
                    upper_flagged_child_coord = max(upper_flagged_child_coord, child_coord)
                end
            end
        end

        if lower_flagged_child_coord > upper_flagged_child_coord # No points flagged or invalid range
            field.regrid_indices = (1, 0) # Or some other conventional "no flagging"
        else
            field.regrid_indices = (lower_flagged_child_coord, upper_flagged_child_coord)
        end
    end

    return nothing
end

"""
STUB: tag_regrid_regions_by_gradient!(fields::Vector{AMRField}, grid::AMRLevel)
Flag regions for regridding based on finite differences of grid functions.
Updates `regrid_indices` in `fields`.
"""
function tag_regrid_regions_by_gradient!(grid::AMRLevel)
    (; ctx, num_grid_points, is_physical_boundary) = grid
    (; fields, buffer_coord, grid_functions_storage_indices) = ctx

    for (field_idx, field) in enumerate(fields)
        if field.pde_type != HYPERBOLIC
            field.regrid_indices = (1, 0) # No flagging
            continue
        end

        storage_idx = grid_functions_storage_indices[field_idx]
        gf = grid.grid_functions_storage[storage_idx]

        lower_flagged_coord = num_grid_points + 1
        upper_flagged_coord = 0

        # Determine effective iteration range on the current `grid`
        iter_start = 1
        # If left boundary is not physical (i.e., interior), apply buffer
        if !is_physical_boundary[1]
            iter_start += buffer_coord
        end

        iter_end = num_grid_points
        # If right boundary is not physical (i.e., interior), apply buffer
        if !is_physical_boundary[2]
            iter_end -= buffer_coord
        end

        # Loop from iter_start to iter_end-1 to compare gf[j] and gf[j+1]
        # So, j can go up to iter_end-1, ensuring j+1 (which is iter_end) is valid.
        # The C loop `for (int jC=start_jC; jC<end_jC-1; jC++)` means jC stops at end_jC-2.
        # If end_jC is the C equivalent of num_grid_points (0-indexed N-1), then end_jC-1 is N-2.
        # Then jC < N-2 means jC goes up to N-3. Accesses gf[N-3] and gf[N-2].
        # Let's adjust Julia: loop j from iter_start to (iter_end - 1)
        if iter_start >= iter_end # Not enough points after applying buffer
            field.regrid_indices = (1, 0)
            continue
        end

        for j_child in iter_start:(iter_end - 1)
            # Ensure j_child and j_child+1 are valid indices for gf
            if j_child < 1 || (j_child + 1) > num_grid_points
                # This check should ideally not be needed if iter_start/iter_end are correct
                continue
            end

            trunc_err = abs(gf[j_child + 1] - gf[j_child])

            if trunc_err > grid.ctx.trunc_err_tolerance
                if lower_flagged_coord > num_grid_points # First time flagging
                    lower_flagged_coord = j_child
                    upper_flagged_coord = j_child
                else
                    lower_flagged_coord = min(lower_flagged_coord, j_child)
                    upper_flagged_coord = max(upper_flagged_coord, j_child)
                end
            end
        end

        if lower_flagged_coord > upper_flagged_coord
            field.regrid_indices = (1, 0)
        else
            # The flagged region is inclusive of [lower_flagged_coord, upper_flagged_coord]
            # If we flag based on cell j (diff between j and j+1), the region needing refinement
            # might include both j and j+1. The C code flags jC.
            # Let's stick to flagging j_child as the start of the difference.
            field.regrid_indices = (lower_flagged_coord, upper_flagged_coord)
        end
    end
    return nothing
end

"""
STUB: determine_overall_regrid_coords!(grid::AMRLevel)
Determine the overall coordinates for creating a new finer grid based on flagged regions
from all fields. Updates `grid.regrid_indices`.
"""
function determine_overall_regrid_coords!(grid::AMRLevel)
    (; ctx, num_grid_points, is_physical_boundary) = grid
    (; fields, buffer_coord) = ctx

    # Initialize with values that will be overridden by any valid flagged region
    # These are 1-based indices for the current `grid`
    overall_lower_coord_jl = num_grid_points + 1
    overall_upper_coord_jl = 0

    # Find min and max coords from all hyperbolic fields
    for field in fields
        if field.pde_type == HYPERBOLIC
            field_lower, field_upper = field.regrid_indices
            # Check if the field has a valid flagged region (lower <= upper and lower >= 1)
            if field_lower >= 1 &&
                field_lower <= field_upper &&
                field_upper <= num_grid_points
                overall_lower_coord_jl = min(overall_lower_coord_jl, field_lower)
                overall_upper_coord_jl = max(overall_upper_coord_jl, field_upper)
            end
        end
    end

    # If no valid region was flagged by any hyperbolic field
    if overall_lower_coord_jl > overall_upper_coord_jl
        grid.regrid_indices = (1, 0) # Indicate no regridding needed
        return nothing
    end

    # Convert to 0-indexed for easier comparison with C logic, then convert back
    lower_c = overall_lower_coord_jl - 1
    upper_c = overall_upper_coord_jl - 1
    nx_c = num_grid_points # C's Nx is number of points, so Nx-1 is max index

    # Apply buffer logic from C code (determine_grid_coords)
    # Note: C variables: lower_child_grid_coord, upper_child_grid_coord, Nx, grid->perim_interior

    # Condition 1: Physical left boundary
    if is_physical_boundary[1] # C: grid->perim_interior[0] == false
        if lower_c < buffer_coord
            lower_c = 0
        end
    end

    # Condition 2: Physical right boundary
    if is_physical_boundary[2] # C: grid->perim_interior[1] == false
        if upper_c > (nx_c - 1 - buffer_coord)
            upper_c = nx_c - 1
        end
    end

    # Condition 3: Interior left boundary
    # C: if ((grid->perim_interior[0]==true) && (lower_child_grid_coord - buffer_coord > 0))
    if !is_physical_boundary[1] # C: grid->perim_interior[0] == true
        if (lower_c - buffer_coord) > 0
            lower_c -= buffer_coord
            # else: C code does not change lower_c if (lower_c - buffer_coord) <= 0
        end
    end

    # Condition 4: Interior right boundary
    # C: if ((grid->perim_interior[1]==true) && (upper_child_grid_coord + buffer_coord < Nx-1))
    # Assuming typo in C code fixed from lower_child_grid_coord to upper_child_grid_coord for this condition.
    if !is_physical_boundary[2] # C: grid->perim_interior[1] == true
        if (upper_c + buffer_coord) < (nx_c - 1)
            upper_c += buffer_coord
            # else: C code does not change upper_c if (upper_c + buffer_coord) >= Nx-1
        end
    end

    # Convert back to 1-based Julia indices and store
    final_lower_jl = lower_c + 1
    final_upper_jl = upper_c + 1

    # Ensure final coordinates are valid and within grid bounds, though C logic should ensure this.
    final_lower_jl = max(1, min(final_lower_jl, num_grid_points))
    final_upper_jl = max(1, min(final_upper_jl, num_grid_points))

    if final_lower_jl > final_upper_jl
        grid.regrid_indices = (1, 0) # Should not happen if initial overall_lower <= overall_upper
    else
        grid.regrid_indices = (final_lower_jl, final_upper_jl)
    end

    return nothing
end

"""
STUB: regrid_level!(grid::AMRLevel)
Orchestrates the regridding process for the level above `grid`.
This involves flagging, determining new grid coordinates, creating the new grid,
interpolating, and injecting old data if applicable.
"""
function regrid_level!(current_L_grid::AMRLevel)
    fields = current_L_grid.ctx.fields
    min_grid_size = current_L_grid.ctx.min_grid_size

    # 1. Flagging based on chosen method
    # If current_L_grid has a parent, we can use Richardson error with parent and current_L_grid.
    # The flags produced will be on current_L_grid.
    if current_L_grid.parent !== nothing
        tag_regridding_regions_by_richardson!(current_L_grid.parent, current_L_grid)
    else
        # Base grid or no parent available, use difference on current_L_grid itself.
        tag_regrid_regions_by_gradient!(current_L_grid)
    end

    # 2. Determine overall regrid coordinates on current_L_grid based on all field flags.
    # Result is stored in current_L_grid.regrid_indices.
    determine_overall_regrid_coords!(current_L_grid)

    # These are the proposed start/end indices on current_L_grid (as parent) for the new child.
    new_child_lower_on_parent, new_child_upper_on_parent = current_L_grid.regrid_indices

    # 3. Validation of new child coordinates & Early Exit if not viable
    if new_child_lower_on_parent > new_child_upper_on_parent ||
        (new_child_upper_on_parent - new_child_lower_on_parent + 1) < min_grid_size

        # If proposed region is invalid/too small, remove existing child if it has no children.
        if current_L_grid.child !== nothing && current_L_grid.child.child === nothing
            amr_remove_level!(current_L_grid.child)
            current_L_grid.child = nothing # Ensure parent's child link is cleared
        end
        # @warn "Regridding aborted: Proposed child grid region is invalid or too small."
        return nothing
    end

    # 4. Handle existing child (old_child)
    old_child = current_L_grid.child

    # 5. Compare with old_child's coordinates to avoid minimal changes
    if old_child !== nothing
        old_child_lower_on_parent = old_child.parent_indices[1]
        old_child_upper_on_parent = old_child.parent_indices[2]
        if abs(old_child_lower_on_parent - new_child_lower_on_parent) <
           current_L_grid.ctx.min_shift_distance &&
            abs(old_child_upper_on_parent - new_child_upper_on_parent) <
           current_L_grid.ctx.min_shift_distance
            # @warn "Regridding aborted: New child grid is too similar to the existing child grid."
            return nothing # New grid is too similar to the old one
        end
    end

    # 6. TODO: Implement C's `make_coords_fit_into_parent_and_fit_grandchild`.
    # This function in C adjusts `new_child_lower_on_parent`, `new_child_upper_on_parent`.
    # For now, we proceed with the coordinates from `determine_overall_regrid_coords!`.
    # After this (if implemented), a re-check of min_grid_size might be needed.
    # if (new_child_upper_on_parent - new_child_lower_on_parent + 1) < min_grid_size:
    #    ... (handle removal of child if it exists and has no children) ...
    #    return nothing

    # 7. Create the new child AMRLevel
    # current_L_grid acts as the parent for this new level.
    new_child = AMRLevel(
        current_L_grid, new_child_lower_on_parent, new_child_upper_on_parent
    )

    # 8. Prolong data from parent (current_L_grid) to new_child
    # Temporarily link new_child to parent for amr_prolong_to_child! to work
    original_actual_child_link = current_L_grid.child # Save current link (could be old_child)
    current_L_grid.child = new_child
    # new_child.parent is already set by its constructor to current_L_grid
    amr_prolong_to_child!(current_L_grid) # Prolongs data into new_child.grid_functions_storage
    current_L_grid.child = original_actual_child_link # Restore parent's original child link for now

    # 9. Handle hierarchy update (replace old_child with new_child)
    grand_child = nothing
    if old_child !== nothing
        @warn "TODO: Implement inject_old_data!(old_child, new_child) to copy overlapping data from old_child to new_child."
        grand_child = old_child.child # Get children of old_child

        # Detach old_child from its parent (current_L_grid) and its child (grand_child)
        # to prepare for its effective removal from this specific parent-child chain.
        old_child.parent = nothing
        old_child.child = nothing
        # old_child itself is not destroyed from memory here, GC will handle if unreferenced.
    end

    # Link new_child into the hierarchy
    new_child.child = grand_child # new_child adopts old_child's children
    if grand_child !== nothing
        grand_child.parent = new_child
        # Update levels of grand_child and its descendants
        current = grand_child
        while current !== nothing
            if current.parent === nothing # Should not happen if logic is correct
                @error "Error in level update: grand_child's parent is nothing."
                break
            end
            current.level = current.parent.level + 1
            current = current.child
        end
    end

    current_L_grid.child = new_child # Parent (current_L_grid) now points to new_child
    # new_child.parent was already set to current_L_grid by AMRLevel constructor.

    # 10. Smooth the parent grid (current_L_grid)
    smooth_all_grid_funcs!(fields, current_L_grid)

    return nothing
end
