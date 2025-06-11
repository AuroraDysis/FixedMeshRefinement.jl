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

    """
        OutputHDF5(filepath::String, grid)

    Create an `OutputHDF5` object. This will create a new HDF5 file (or
    truncate an existing one) and set up the necessary groups and extendible
    datasets based on the initial grid structure.
    The datasets for fields (`psi`, `Pi`) will have dimensions `(time, x)`
    which is convenient for reading with Python's h5py library.
    """
    function OutputHDF5(
        filepath::String, grid::Grid{NumState,NumDiagnostic}
    ) where {NumState,NumDiagnostic}
        # Ensure the directory exists
        dir = dirname(filepath)
        if !isdir(dir)
            mkpath(dir)
        end

        (; num_levels, levels) = grid

        local file::HDF5.File
        if isfile(filepath)
            file = h5open(filepath, "r+")
        else
            file = h5open(filepath, "w")

            try
                for l in 1:num_levels
                    level = levels[l]

                    x = parent(get_x(level))
                    num_points = length(x)

                    g = create_group(file, "level$(lpad(l, 2, '0'))")

                    # Create static dataset for coordinates
                    write(g, "x", collect(x))

                    # Create extendible datasets for time-dependent data
                    t_ds_space = HDF5.dataspace((3, 0); max_dims=(3, -1))
                    field_ds_space = HDF5.dataspace(
                        (num_points, NumState, 0); max_dims=(num_points, NumState, -1)
                    )
                    diagnostic_ds_space = HDF5.dataspace(
                        (num_points, NumDiagnostic, 0);
                        max_dims=(num_points, NumDiagnostic, -1),
                    )

                    # Use chunking for performance with extendible datasets
                    t_chunk = (3, 1024)
                    field_chunk = (num_points, NumState, 1) # chunked by time slice
                    diagnostic_chunk = (num_points, NumDiagnostic, 1)

                    create_dataset(g, "t", Float64, t_ds_space; chunk=t_chunk)
                    create_dataset(g, "state", Float64, field_ds_space; chunk=field_chunk)
                    create_dataset(
                        g,
                        "diagnostic",
                        Float64,
                        diagnostic_ds_space;
                        chunk=diagnostic_chunk,
                    )
                end
            catch e
                close(file)
                rethrow(e)
            end
        end

        level_grid_points = get_maximum_grid_points.(levels)
        return new(filepath, file, level_grid_points)
    end
end

"""
    append_data(out::OutputHDF5, grid)

Append a new time slice of data from the grid to the HDF5 file.
"""
function append_data(
    out::OutputHDF5, grid::Grid{NumState,NumDiagnostic}
) where {NumState,NumDiagnostic}
    (; num_levels, levels) = grid

    for l in 1:num_levels
        level = levels[l]
        num_points = out.level_grid_points[l]
        state = get_state(level)
        diagnostic = get_diagnostic_state(level)

        # map indices to the HDF5 dataset
        interior_indices = get_interior_indices(level)
        offset_indices = get_offset_indices(level)
        nleft = first(interior_indices) - first(offset_indices)
        dset_interior_indices = interior_indices .+ nleft
        dset_offset_indices = offset_indices .+ nleft
        dset_left_indices = first(dset_offset_indices):(first(dset_interior_indices) - 1)
        dset_right_indices = (last(dset_interior_indices) + 1):last(dset_offset_indices)

        g = out.file["level$(lpad(l, 2, '0'))"]

        # Append time
        dset_t = g["t"]
        current_len = size(dset_t, 2)
        new_len = current_len + 1
        HDF5.set_extent_dims(dset_t, (3, new_len))
        dset_t[1, new_len] = level.t
        dset_t[2, new_len] = first(dset_interior_indices) # first interior index
        dset_t[3, new_len] = last(dset_interior_indices) # last interior index

        # Append state
        dset_state = g["state"]
        HDF5.set_extent_dims(dset_state, (num_points, NumState, new_len))
        dset_state[dset_left_indices, :, new_len] .= NaN
        dset_state[dset_interior_indices, :, new_len] .= @view(state[interior_indices, :])
        dset_state[dset_right_indices, :, new_len] .= NaN

        # Append diagnostic
        dset_diagnostic = g["diagnostic"]
        HDF5.set_extent_dims(dset_diagnostic, (num_points, NumDiagnostic, new_len))
        dset_diagnostic[dset_left_indices, :, new_len] .= NaN
        dset_diagnostic[dset_interior_indices, :, new_len] .= @view(
            diagnostic[interior_indices, :]
        )
        dset_diagnostic[dset_right_indices, :, new_len] .= NaN
    end
end

"""
    Base.close(out::OutputHDF5)

Close the HDF5 file associated with the `OutputHDF5` object.
"""
function Base.close(out::OutputHDF5)
    return close(out.file)
end
