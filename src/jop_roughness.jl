"""
    A = JopRoughness(sp, dim, w)

`A`  can be used to regularize an optimization problem, applying a penalty to
models that are non-smooth.  This is similar to a finite difference operator,
but here no care is taken to ensure that `A` computes a derivative.  `dim` is
the dimension that one wishes to smooth and `w` is the half-width of the smoother.
The width of the window determines how strong the penalty is for being non-smooth.

For example, the form of `A` for `w=1` is,
```julia
A=
[
(1/3-2/3) 1/3       0         0         ;
1/3       (1/3-3/3) 1/3       0         ;
0         1/3       (1/3-3/3) 1/3       ;
0         0         1/3       (1/3-2/3)
]
```
"""
JopRoughness(sp::JetSpace, dim::Int; w::Int=1) = JopLn(dom = sp, rng = sp, df! = JopRoughness_df!, s = (dim=dim, w=w))

export JopRoughness

function JopRoughness_df!(d, m; dim, w, kwargs...)
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)
    roughness_helper(d, m, w, n, Rpre, Rpost)
end

function roughness_helper(d::AbstractArray{T}, m::AbstractArray{T}, w, n, Rpre, Rpost) where {T}
    d .= 0
    _w = 2*w + 1
    for Ipost in Rpost
        for i = 1:n
            f = -one(T)
            for j = 1:_w
                k = i + j - 1 - w
                if k ∈ 1:n
                    for Ipre in Rpre
                        d[Ipre, i, Ipost] += m[Ipre, k, Ipost]
                    end
                else
                    f += one(T)/_w
                end
            end
            for Ipre in Rpre
                d[Ipre, i, Ipost] /= _w
                d[Ipre, i, Ipost] += m[Ipre, i, Ipost]*f
            end
        end
    end
    d
end
