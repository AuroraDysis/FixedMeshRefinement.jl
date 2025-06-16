using HDF5

export OutputHDF5, append_data

"""
    OutputHDF5

A struct to handle writing data to an HDF5 file slice by slice.
This is designed to create a single HDF5 file for a simulation run,
with datasets that are extended over time.
"""
mutable struct OutputHDF5
    path::String
    file::HDF5.File
    level_grid_points::Vector{Int}
    save_merged_data::Bool
    max_merged_length::Int
    merged_vec::Vector{Float64}

    """
        OutputHDF5(filepath::String, grid)

    Create an `OutputHDF5` object. This will create a new HDF5 file (or
    truncate an existing one) and set up the necessary groups and extendible
    datasets based on the initial grid structure.
    """
    function OutputHDF5(filepath::String, grid::Grid; save_merged_data::Bool=false)
        # Ensure the directory exists
        dir = dirname(filepath)
        if !isdir(dir)
            mkpath(dir)
        end

        (; num_state_variables, num_diagnostic_variables) = grid
        num_levels = get_num_levels(grid)
        max_merged_length = -1

        local file::HDF5.File
        if isfile(filepath)
            file = h5open(filepath, "r+")
        else
            file = h5open(filepath, "w")

            try
                for l in 1:num_levels
                    level = get_level(grid, l)

                    x = parent(get_x(level))
                    num_points = length(x)

                    g = create_group(file, "level$(lpad(l, 2, '0'))")

                    # Create static dataset for coordinates
                    write(g, "x", collect(x))

                    # Create extendible datasets for time-dependent data
                    t_ds_space = HDF5.dataspace((3, 0); max_dims=(3, -1))
                    t_chunk = (3, 1024)
                    create_dataset(g, "t", Float64, t_ds_space; chunk=t_chunk)

                    state_ds_space = HDF5.dataspace(
                        (num_points, num_state_variables, 0);
                        max_dims=(num_points, num_state_variables, -1),
                    )
                    state_chunk = (num_points, num_state_variables, 1) # chunked by time slice
                    create_dataset(g, "state", Float64, state_ds_space; chunk=state_chunk)

                    if num_diagnostic_variables > 0
                        diagnostic_ds_space = HDF5.dataspace(
                            (num_points, num_diagnostic_variables, 0);
                            max_dims=(num_points, num_diagnostic_variables, -1),
                        )
                        diagnostic_chunk = (num_points, num_diagnostic_variables, 1)
                        create_dataset(
                            g,
                            "diagnostic",
                            Float64,
                            diagnostic_ds_space;
                            chunk=diagnostic_chunk,
                        )
                    end
                end

                if save_merged_data
                    g_merged = create_group(file, "merged_state")

                    x, y = merge_grid_levels(grid, l -> get_state(l)[:, 1])
                    write(g_merged, "x", collect(x))

                    t_merged_ds_space = HDF5.dataspace((0,); max_dims=(-1,))
                    t_merged_chunk = (1024,)
                    create_dataset(g_merged, "t", Float64, t_merged_ds_space; chunk=t_merged_chunk)

                    max_merged_length = max(max_merged_length, length(x))
                    merged_ds_space = HDF5.dataspace((max_merged_length, 0); max_dims=(max_merged_length, -1))
                    merged_chunk = (max_merged_length, 1)
                    create_dataset(g_merged, "psi", Float64, merged_ds_space; chunk=merged_chunk)
                end
            catch e
                close(file)
                rethrow(e)
            end
        end

        level_grid_points = zeros(Int, num_levels)
        for l in 1:num_levels
            level = get_level(grid, l)
            level_grid_points[l] = get_maximum_grid_points(level)
        end
        merged_vec = zeros(max_merged_length)
        return new(filepath, file, level_grid_points, save_merged_data, max_merged_length, merged_vec)
    end
end

"""
    append_data(out::OutputHDF5, grid)

Append a new time slice of data from the grid to the HDF5 file.
"""
function append_data(out::OutputHDF5, grid::Grid)
    (; num_state_variables, num_diagnostic_variables) = grid
    num_levels = get_num_levels(grid)

    for l in 1:num_levels
        level = get_level(grid, l)
        num_points = out.level_grid_points[l]
        state = get_state(level)
        diagnostic = get_diagnostic_state(level)

        # map indices to the HDF5 dataset
        interior_indices = get_interior_indices(level)
        offset_indices = get_offset_indices(level)
        left_indices = first(offset_indices):(first(interior_indices) - 1)
        right_indices = (last(interior_indices) + 1):last(offset_indices)
        state[left_indices, :] .= NaN
        state[right_indices, :] .= NaN

        # interior indices for the HDF5 dataset
        nleft = first(interior_indices) - first(offset_indices)
        dset_interior_indices = interior_indices .+ nleft

        g = out.file["level$(lpad(l, 2, '0'))"]

        # Append time and interior indices
        dset_t = g["t"]
        current_len = size(dset_t, 2)
        new_len = current_len + 1
        HDF5.set_extent_dims(dset_t, (3, new_len))
        dset_t[1, new_len] = level.t
        dset_t[2, new_len] = first(dset_interior_indices) # first interior index
        dset_t[3, new_len] = last(dset_interior_indices) # last interior index

        # Append state
        dset_state = g["state"]
        HDF5.set_extent_dims(dset_state, (num_points, num_state_variables, new_len))
        dset_state[:, :, new_len] = parent(state)

        # Append diagnostic
        dset_diagnostic = g["diagnostic"]
        HDF5.set_extent_dims(
            dset_diagnostic, (num_points, num_diagnostic_variables, new_len)
        )
        dset_diagnostic[:, :, new_len] = parent(diagnostic)
    end

    if out.save_merged_data
        g_merged = out.file["merged_state"]

        dset_merged_t = g_merged["t"]
        current_len = size(dset_merged_t, 1)
        new_len = current_len + 1
        HDF5.set_extent_dims(dset_merged_t, (1, new_len))
        dset_merged_t[new_len] = grid.t

        dset_merged = g_merged["psi"]
        HDF5.set_extent_dims(dset_merged, (out.max_merged_length, new_len))
        x, y = merge_grid_levels(grid, l -> get_state(l)[:, 1])
        lo_idx = out.max_merged_length - length(y) + 1
        out.merged_vec[1:lo_idx-1] .= NaN
        out.merged_vec[lo_idx:end] .= y
        dset_merged[:, new_len] = out.merged_vec
    end
end

"""
    Base.close(out::OutputHDF5)

Close the HDF5 file associated with the `OutputHDF5` object.
"""
function Base.close(out::OutputHDF5)
    return close(out.file)
end
