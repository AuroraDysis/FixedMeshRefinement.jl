export typetol, extrapolate!

"""
    typetol(TF::Type{<:AbstractFloat})

Return a tolerance value for a given float type `TF`. The tolerance is calculated as
`eps(TF)^(4/5)`. This is useful for setting default tolerances in numerical algorithms
where machine precision is a factor.

# Arguments
- `TF::Type{<:AbstractFloat}`: The floating point type.

# Returns
- A tolerance value of the same floating point type.

# Examples
```jldoctest
julia> typetol(Float64)
3.666852862501036e-11

julia> typetol(Float32)
2.422181f-5
```
"""
function typetol(TF::Type{<:AbstractFloat})
    return eps(TF)^(2//3)
end

function isapprox_tol(a::T, b::T) where {T<:AbstractFloat}
    tol = typetol(T)
    return isapprox(a, b; atol=tol, rtol=tol)
end

function extrapolate!(out, input, order::Integer)
    length(input) == order || error("Length of input must be equal to order.")

    if order == 1
        @.. out = input[1]
    elseif order == 2
        @.. out = 2 * input[1] - input[2]
    elseif order == 3
        @.. out = 3 * input[1] - 3 * input[2] + input[3]
    elseif order == 4
        @.. out = 4 * input[1] - 6 * input[2] + 4 * input[3] - input[4]
    elseif order == 5
        @.. out = 5 * input[1] - 10 * input[2] + 10 * input[3] - 5 * input[4] + input[5]
    elseif order == 6
        @.. out =
            6 * input[1] - 15 * input[2] + 20 * input[3] - 15 * input[4] + 6 * input[5] -
            input[6]
    else
        error("Extrapolation order $order is not supported.")
    end
end
