function typetol(TF::Type{<:AbstractFloat})
    return eps(TF)^(4//5)
end
