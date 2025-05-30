
"""
amr_init_grid_hierarchy(fields_vector; num_total_time_steps, num_base_spatial_points, save_interval, cfl, domain_bounds, use_excision) -> AMRContext
Allocate base grid (level 0) and an initial level-1 grid refined everywhere.
"""
function amr_init_grid_hierarchy(
    fields_vector::Vector{AMRField};
    num_total_time_steps::Int,
    num_base_spatial_points::Int,
    save_interval::Int,
    cfl::Float64,
    domain_bounds::Tuple{Float64,Float64},
    use_excision::Bool,
)
    # Convert linked list of fields to a Vector and compute total num_gfs
    total_grid_functions = sum(
        field.num_time_levels + field.num_extrapolation_levels for field in fields_vector
    )

    # Base (shadow) grid -------------------------------------------------------
    base_grid_obj = AMRGrid(0, num_base_spatial_points, total_grid_functions)
    base_grid_obj.bounding_box = MVector(domain_bounds...)
    base_grid_obj.dx = (domain_bounds[2] - domain_bounds[1]) / (num_base_spatial_points - 1)
    base_grid_obj.dt = cfl * base_grid_obj.dx
    base_grid_obj.parent_grid_bounds = MVector(0, num_base_spatial_points - 1)
    base_grid_obj.is_physical_boundary = MVector(false, false) # Assuming physical boundaries are at the ends of the domain
    base_grid_obj.excision_on = use_excision

    # Create ctx first
    ctx = AMRContext(
        cfl, num_total_time_steps, save_interval, fields_vector, base_grid_obj, use_excision
    )

    # Add level-1 grid refined everywhere initially
    level_one_grid = amr_make_finer_grid!(
        ctx, base_grid_obj, 0, num_base_spatial_points - 1
    )
    base_grid_obj.child = level_one_grid

    return ctx
end