using HDF5
using Printf

export OutputHDF5, append_data

"""
    OutputHDF5(output_dir::String, file_prefix::String, grid::Grid)

A struct to handle writing data to separate HDF5 files for each time step.
This creates multiple small HDF5 files instead of one large file.
Also maintains a CSV index table of simulation snapshots for easy reference.
"""
mutable struct OutputHDF5
    output_dir::String
    file_prefix::String
    level_grid_points::Vector{Int}
    num_state_variables::Int
    num_diagnostic_variables::Int
    snapshot_index::OutputCSV

    """
        OutputHDF5(output_dir::String, file_prefix::String, grid)

    Create an `OutputHDF5` object. This will create the output directory
    and prepare for writing separate HDF5 files for each time step.
    Also creates a CSV index table to track simulation snapshots.
    """
    function OutputHDF5(output_dir::String, file_prefix::String, grid::Grid)
        # Ensure the output directory exists
        if !isdir(output_dir)
            mkpath(output_dir)
        end

        (; num_state_variables, num_diagnostic_variables) = grid
        num_levels = get_num_levels(grid)

        level_grid_points = zeros(Int, num_levels)
        for l in 1:num_levels
            level = get_level(grid, l)
            level_grid_points[l] = get_maximum_grid_points(level)
        end

        # Create CSV index table for snapshot tracking
        csv_path = joinpath(output_dir, "_index.csv")
        snapshot_index = OutputCSV(csv_path, ["step", "time", "filename"])

        # Create the OutputHDF5 object
        return new(
            output_dir,
            file_prefix,
            level_grid_points,
            num_state_variables,
            num_diagnostic_variables,
            snapshot_index,
        )
    end
end

"""
    append_data(out::OutputHDF5, grid)

Save the current state of the grid to a new HDF5 file.
Each call creates a separate file with the current time step data.
Also records the snapshot information in the CSV index table.
"""
function append_data(out::OutputHDF5, grid::Grid, iter)
    # Create filename with counter and time
    t = grid.t
    filename = "$(out.file_prefix)_$(lpad(iter, 8, '0'))_t_$(Printf.@sprintf("%.3f", t)).h5"
    full_path = joinpath(out.output_dir, filename)

    num_levels = get_num_levels(grid)

    # Create new file for this time step
    file = h5open(full_path, "w")

    try
        # Write global attributes
        HDF5.attributes(file)["time"] = t
        HDF5.attributes(file)["step"] = iter
        HDF5.attributes(file)["num_levels"] = num_levels

        for l in 1:num_levels
            level = get_level(grid, l)
            state = get_state(level)
            diagnostic = get_diagnostic_state(level)

            # map indices to the HDF5 dataset
            interior_indices = get_interior_indices(level)
            offset_indices = get_offset_indices(level)
            left_indices = first(offset_indices):(first(interior_indices) - 1)
            right_indices = (last(interior_indices) + 1):last(offset_indices)

            state[:, left_indices] .= NaN
            state[:, right_indices] .= NaN
            diagnostic[:, left_indices] .= NaN
            diagnostic[:, right_indices] .= NaN

            # interior indices for the HDF5 dataset
            nleft = first(interior_indices) - first(offset_indices)
            dset_interior_indices = interior_indices .+ nleft

            # Create group for this level
            g = create_group(file, "level$(lpad(l, 2, '0'))")

            # Write level attributes
            HDF5.attributes(g)["level_number"] = l
            HDF5.attributes(g)["time"] = level.t
            HDF5.attributes(g)["first_interior_index"] = first(dset_interior_indices)
            HDF5.attributes(g)["last_interior_index"] = last(dset_interior_indices)

            # Write coordinate data
            x = parent(get_x(level))
            write(g, "x", collect(x))

            # Write state data
            write(g, "state", parent(state))

            # Write diagnostic data
            write(g, "diagnostic", parent(diagnostic))
        end
    finally
        close(file)
    end

    # Record snapshot information in index table
    write_row(out.snapshot_index, (iter, t, filename))

    return full_path
end

"""
    Base.close(out::OutputHDF5)

Close the CSV snapshot index associated with the OutputHDF5 object.
"""
function Base.close(out::OutputHDF5)
    return close(out.snapshot_index)
end
