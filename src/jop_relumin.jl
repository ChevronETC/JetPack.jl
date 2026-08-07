"""
    F = JopReluMin(spc)

where `F` is the relu negative activation function operator with domain and range given by `spc::JetSpace`.

Note that `F` is not a liner function in the sense that for general vectors `x` and `y` 
```
    F(α x + y) != α F(x) + F(y)
```
However if `x` and `y` are either both strictly positive, or both strictly negative, F is linear.
"""
JopReluMin(spc::JetSpace) = JopLn(dom = spc, rng = spc, df! = JopReluMin_df!)
export JopReluMin

JopReluMin_df!(d::AbstractArray{T}, m::AbstractArray{T}; kwargs...) where {T <: Real} = d .= map(x -> x < 0 ? x : 0, m)
