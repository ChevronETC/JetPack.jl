"""
    A = JopRestriction(dom[, rng], indices)

Apply a restriction operator mapping from `dom::JetSpace` to
`rng::JetSpace`, and using a one-dimensional array of indices
where `indices` is one of:

* `indices::Vector{Int,1}` - for restriction of 1-D arrays
* `indices::Vector{NTuple{N,Int}}` - for restriction of N-D arrays (N>1)

Note that in the case where the domain is N-dimension, the range is always stored
using a 1-dimension array.

# example (1-D):
```julia
A = JopRestriction(JotSpace(Float64,4),[1,3])
m = [1.0;2.0;3.0;4.0]
d = A*m # d=[1.0;3.0]
```

# example (2-D):
```julia
A = JopRestriction(JotSpace(Float64,2,2),[(1,2),(2,2)])
m = [1.0 2.0;3.0 4.0]
d = A*m # d=[2.0,4.0]
```
"""
function JopRestriction(dom::JetAbstractSpace{T}, indices::Vector) where {T}
    rng = JetSpace(T, length(indices))
    JopRestriction(dom, rng, indices)
end

function JopRestriction(dom::JetAbstractSpace, rng::JetAbstractSpace, indices::Vector)
    _indices = [CartesianIndex((index...,)) for index in indices]
    JopRestriction(dom, rng, _indices)
end

function JopRestriction(dom::JetAbstractSpace, indices::Vector{C}) where {C<:CartesianIndex}
    rng = JetSpace(T, length(indices))
    JopRestriction(dom, rng, indices)
end

function JopRestriction(dom::JetAbstractSpace, rng::JetAbstractSpace, indices::Vector{C}) where {C<:CartesianIndex}
    JopLn(dom = dom, rng = rng, df! = JopRestriction_df!, df′! = JopRestriction_df′!, s = (indices=indices,))
end

export JopRestriction

function JopRestriction_df!(d, m; indices, kwargs...)
    for i = 1:length(indices)
        d[i] = m[indices[i]]
    end
    d
end

function JopRestriction_df′!(m, d; indices, kwargs...)
    m .= 0
    for i = 1:length(indices)
        m[indices[i]] += d[i]
    end
    m
end
