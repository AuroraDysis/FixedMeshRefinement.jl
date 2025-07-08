export OutputCSV, write_row

"""
    OutputCSV

A struct to handle writing data to a CSV file line by line.
"""
mutable struct OutputCSV
    path::String
    io::IO
    num_columns::Int
    header::Union{Vector{String},Nothing}

    """
        OutputCSV(filepath::String)

    Create an `OutputCSV` object and open the file for writing.
    """
    function OutputCSV(filepath::String, num_columns::Int)
        if isfile(filepath)
            error("File $filepath already exists, please delete it first")
        end

        io = open(filepath, "w")
        return new(filepath, io, num_columns, nothing)
    end

    function OutputCSV(filepath::String, header::Vector{String})
        if isfile(filepath)
            error("File $filepath already exists, please delete it first")
        end

        io = open(filepath, "w")
        println(io, join(header, ","))
        return new(filepath, io, length(header), header)
    end
end

"""
    write_row(out::OutputCSV, row)

Write a single row to the CSV file. The `row` must be a `Tuple` or `Vector`.
The header is automatically written on the first call based on the names in the `row`.
"""
function write_row(out::OutputCSV, row)
    return println(out.io, join(row, ","))
end

"""
    Base.close(out::OutputCSV)

Close the file stream associated with the `OutputCSV` object. It's safe to call this
multiple times.
"""
function Base.close(out::OutputCSV)
    if isopen(out.io)
        close(out.io)
    end
end
