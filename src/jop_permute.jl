"""
    A = JopPermute(sp, dims, perm)

where `sp::JetSpace` is the domain and range of the operator, `dims::NTuple` is a tuple of dimensions
that one wishes to permute, and `perm` is a tuple of arrays with permutation indices.

# Example 1:
```julia
A = JopPermute(JetSpace(Float64,4,2),(1,),([3;1;2;4],))
m = [11 12 ; 21 22 ; 31 32 ; 41 42]
d = A*m # d=[31 32 ; 11 12 ; 21 22 ; 41 42]
```

# Example 2:
```julia
A = JopPermute(JetSpace(Float64,3,2),(1,2),([3;2;1],[2;1]))
m = [11 12 ; 21 22 ; 31 32]
d = A*m # d = [32 31 ; 22 21 ; 12 11]
```

# Notes:
Currently, this is only implmented for 1D, 2D and 3D arrays.
"""
JopPermute(sp, dims, perm) = JopLn(dom = sp, rng = sp, df! = JopPermute_df!, s = (dims=dims, perm=perm))

export JopPermute

# 1D
function JopPermute_df!(d::AbstractVector, m::AbstractVector; perm, kwargs...)
    d[:] = m[perm[1]]
    d
end

# 2D
function JopPermute_df!(d::AbstractMatrix, m::AbstractMatrix; perm, dims, kwargs...)
    d .= m
    for i=1:length(dims)
        if dims[i] == 1
            d .= d[perm[i],:]
        end
        if dims[i] == 2
            d .= d[:,perm[i]]
        end
    end
    d
end

# 3D
function JopPermute_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; perm, dims, kwargs...) where {T}
    d .= m
    for i=1:length(dims)
        if dims[i] == 1
            d .= d[perm[i],:,:]
        end
        if dims[i] == 2
            d .= d[:,perm[i],:]
        end
        if dims[i] == 3
            d .= d[:,:,perm[i]]
        end
    end
    d
end
