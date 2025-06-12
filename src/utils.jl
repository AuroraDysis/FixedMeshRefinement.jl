export typetol

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
