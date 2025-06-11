module CsvOutput

using StaticArrays

export OutputCSV, write_row

"""
    OutputCSV

A struct to handle writing data to a CSV file line by line.
"""
mutable struct OutputCSV{NumColumns}
    path::String
    io::IO
    header::Union{SVector{NumColumns,String},Nothing}

    """
        OutputCSV(filepath::String)

    Create an `OutputCSV` object and open the file for writing.
    """
    function OutputCSV(filepath::String)
        io = open(filepath, "a")
        return new{0}(filepath, io, nothing)
    end

    function OutputCSV(filepath::String, header::SVector{NumColumns,String}) where {NumColumns}
        io = open(filepath, "a")
        println(io, join(header, ","))
        return new{NumColumns}(filepath, io, header)
    end
end

"""
    write_row(out::OutputCSV, row::NamedTuple)

Write a single row to the CSV file. The `row` must be a `NamedTuple`.
The header is automatically written on the first call based on the names in the `row`.
"""
function write_row(
    out::OutputCSV{NumColumns}, row::SVector{NumColumns,Float64}
) where {NumColumns}
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

end # module CsvOutput
