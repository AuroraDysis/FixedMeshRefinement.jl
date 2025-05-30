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
        return new(
            max_levels,
            refinement_ratio,
            regrid_cadence,
            buffer_coord,
            min_grid_size,
            fields,
            base_grid,
            use_excision,
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
    regrid_bounds::NTuple{2,Int}         # (lower, upper) flagged coordinates for regridding

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
    parent_grid_bounds::NTuple{2,Int}                                  # coordinates inside parent grid (inclusive)
    regrid_bounds::NTuple{2,Int}                                       # flagged coordinates for regridding
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
    parent_grid_bounds = (0, 0)
    regrid_bounds = (0, 0)
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
        parent_grid_bounds,
        regrid_bounds,
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
    parent_grid_bounds_val = (start_idx, end_idx)

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
        parent_grid_bounds_val,             # parent_grid_bounds
        (0, 0),                             # regrid_bounds (default placeholder)
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
amr_restrict_to_parent!(ctx::AMRContext, current_grid::AMRLevel)
Restrict overlapping grid-functions from `current_grid` into its `parent`.
This version assumes the 'index' corresponded to the first data slot of each field.
"""
function amr_restrict_to_parent!(ctx::AMRContext, current_grid::AMRLevel)
    parent = current_grid.parent

    calculated_field_data_slot_index = 0 # 0-based index for the start of current field's data block in grid_functions
    for field_obj in ctx.fields
        # The Julia index (1-based) for the specific Vector{Float64} in grid_functions
        # This assumes the original fld.index pointed to the first time_level (n=0) data for that field.
        actual_julia_index_in_grid_funcs = calculated_field_data_slot_index + 1

        copy_to_parent_grid!(
            parent.grid_functions[actual_julia_index_in_grid_funcs],
            current_grid.grid_functions[actual_julia_index_in_grid_funcs],
            current_grid.parent_grid_bounds[1],
            ctx,
        )

        # Advance the base slot index to the start of the next field's data block
        calculated_field_data_slot_index +=
            field_obj.num_time_levels + field_obj.num_extrapolation_levels
    end
    return nothing
end
