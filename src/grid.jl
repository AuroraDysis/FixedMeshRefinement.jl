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
    base_grid::AMRLevel

    cfl::TF # CFL number

    num_grid_functions::Int
    grid_functions_storage_indices::Vector{Int}

    function AMRContext(
        fields::Vector{AMRField},
        base_grid::AMRLevel,
        use_excision::Bool;
        max_levels::Int=8,
        refinement_ratio::Int=2,
        regrid_cadence::Int=80,
        buffer_coord::Int=40,
        min_grid_size::Int=40,
    )
        num_grid_functions = [
            field.num_time_levels + field.num_extrapolation_levels for field in fields
        ]
        grid_functions_storage_indices = Vector{Int}(undef, num_grid_functions)
        current_index = 1
        for (i, field) in enumerate(fields)
            grid_functions_storage_indices[i] = current_index
            current_index += field.num_time_levels + field.num_extrapolation_levels
        end
        return new(
            max_levels,
            refinement_ratio,
            regrid_cadence,
            buffer_coord,
            min_grid_size,
            fields,
            base_grid,
            use_excision,
            num_grid_functions,
            grid_functions_storage_indices,
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
        (0, 0),                             # regrid_indices
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
        for j in 1:num_points-1
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
    parent_grid = grid_to_remove.parent
    child_grid = grid_to_remove.child

    if !isnothing(parent_grid)
        parent_grid.child = child_grid # Link parent to grandchild
    end

    if !isnothing(child_grid)
        child_grid.parent = parent_grid # Link grandchild to grandparent
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
    # Placeholder for Kreiss-Oliger smoothing logic
    @warn "smooth_all_grid_funcs! is a stub and not implemented."
    return nothing
end

"""
STUB: flag_regridding_regions_richardson!(fields::Vector{AMRField}, parent_grid::AMRLevel, grid::AMRLevel)
Flag regions for regridding based on Richardson extrapolation error.
Updates `regrid_indices` in `fields`.
"""
function flag_regridding_regions_richardson!(
    fields::Vector{AMRField}, parent_grid::AMRLevel, grid::AMRLevel
)
    # Placeholder for Richardson-based flagging
    @warn "flag_regridding_regions_richardson! is a stub and not implemented."
    for field in fields
        field.regrid_indices = (0, 0) # Default no flagging
    end
    return nothing
end

"""
STUB: flag_regridding_regions_difference!(fields::Vector{AMRField}, grid::AMRLevel)
Flag regions for regridding based on finite differences of grid functions.
Updates `regrid_indices` in `fields`.
"""
function flag_regridding_regions_difference!(fields::Vector{AMRField}, grid::AMRLevel)
    # Placeholder for difference-based flagging
    @warn "flag_regridding_regions_difference! is a stub and not implemented."
    for field in fields
        field.regrid_indices = (0, 0) # Default no flagging
    end
    return nothing
end

"""
STUB: determine_overall_regrid_coords!(fields::Vector{AMRField}, grid::AMRLevel)
Determine the overall coordinates for creating a new finer grid based on flagged regions
from all fields. Updates `grid.regrid_indices`.
"""
function determine_overall_regrid_coords!(fields::Vector{AMRField}, grid::AMRLevel)
    # Placeholder for logic to combine flagged regions
    @warn "determine_overall_regrid_coords! is a stub and not implemented."
    grid.regrid_indices = (0, 0) # Default no regridding
    return nothing
end

"""
STUB: regrid_level!(grid::AMRLevel)
Orchestrates the regridding process for the level above `grid`.
This involves flagging, determining new grid coordinates, creating the new grid,
interpolating, and injecting old data if applicable.
"""
function regrid_level!(grid::AMRLevel)
    # Placeholder for the main regridding orchestration logic
    @warn "regrid_level! is a stub and not implemented."
    # This would call flagging functions, determine_overall_regrid_coords!,
    # then potentially amr_refine_level! or amr_insert_level! with new logic.
    return nothing
end
