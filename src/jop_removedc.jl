"""
    JopReoveDC(s::JetAbstractSpace{T,N})

F(x) = x - zero frequency component of x along the fastest dimension
Uses Fourier transform, as the DC component differs from the mean as a function of interval length and zero padding.
```
"""
function JopRemoveDC(spc::JetAbstractSpace{T,N}) where {T,N}
    JopLn(dom = spc, rng = spc, df! = JopRemoveDC_df!, df′! = JopRemoveDC_df′!)
end

export JopRemoveDC

function JopRemoveDC_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kwargs...) where {T}
    n1 = size(m,1)
    n2 = div(length(m),n1)
    _m = reshape(m,n1,n2)
    _d = reshape(d,n1,n2)
    for k2 = 1:n2
        mean = sum(_m[:,k2]) / length(_m[:,k2])
        _d[:,k2] .= _m[:,k2] .- mean
    end
    _d
end

function JopRemoveDC_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; kwargs...) where {T}
    n1 = size(m,1)
    n2 = div(length(m),n1)
    _m = reshape(m,n1,n2)
    _d = reshape(d,n1,n2)
    for k2 = 1:n2
        mean = sum(_d[:,k2]) / length(_d[:,k2])
        _m[:,k2] .= _d[:,k2] .- mean
    end
    _m
end
