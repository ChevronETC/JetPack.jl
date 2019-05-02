"""
    A = JopGradient(dom, δ)

Gradient for a 2D or 3D arrays. The range will have one more
dimension than the domain and contain the components of the gradient
in each dimension.

* `dom::JotSpace{T,N}` is the domain of the operator.
* `δ::NTuple{T,N}` is the grid spacing in each dimension.

# Examples:

## 2D
```julia
nz,nx = 11,12
dom = JetSpace(Float32, nz, nx)
A = JopGradient(dom, (1.0,1.0))
size(domain(A)) # (11,12)
size(range(A))  # (11,12,2)
```

## 3D
```julia
nz,ny,nx = 11,12,13
dom = JetSpace(Float32, nz, ny, nx)
A = JopGradient(dom,(1.0,1.0,1.0))
size(domain(A)) # (11,12,13)
size(range(A))  # (11,12,13,3)
```
"""
function JopGradient(dom::JetSpace{T,N}, δ::NTuple{N,T}) where {T,N}
    rng = JetSpace(T, (size(dom)...,length(size(dom))))
    if length(size(dom)) != length(δ)
        error("dimensionality of length(size(dom)) ($(length(size(dom)))) and length(δ) ($(length(δ))) do not agree!")
    end
    JopLn(dom = dom, rng = rng, df! = JopGradient_df!, df′! = JopGradient_df′!, s = (δ=δ,))
end

export JopGradient

# ....................................................................................
# 2D
# ....................................................................................
function JopGradient_df!(d::AbstractArray{T,3}, m::AbstractArray{T,2}; δ, kwargs...) where {T}
    d .= 0
    n1,n2 = size(m)
    for k1=1:n1-1, k2=1:n2
        d[k1,k2,1] += m[k1+1,k2] / δ[1]
        d[k1,k2,1] -= m[k1+0,k2] / δ[1]
    end
    for k1=1:n1, k2=1:n2-1
        d[k1,k2,2] += m[k1,k2+1] / δ[2]
        d[k1,k2,2] -= m[k1,k2+0] / δ[2]
    end
    d
end

function JopGradient_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,3}; δ, kwargs...) where {T}
    m .= 0
    n1,n2 = size(m)
    for k1=1:n1-1, k2=1:n2
        m[k1+1,k2] += d[k1,k2,1] / δ[1]
        m[k1+0,k2] -= d[k1,k2,1] / δ[1]
    end
    for k1=1:n1, k2=1:n2-1
        m[k1,k2+1] += d[k1,k2,2] / δ[2]
        m[k1,k2+0] -= d[k1,k2,2] / δ[2]
    end
    m
end

# ....................................................................................
# 3D
# ....................................................................................
function JopGradient_df!(d::AbstractArray{T,4}, m::AbstractArray{T,3}; δ, kwargs...) where {T}
    d .= 0
    n1,n2,n3 = size(d)
    for k1=1:n1-1, k2=1:n2, k3=1:n3
        d[k1,k2,k3,1] += m[k1+1,k2,k3] / δ[1]
        d[k1,k2,k3,1] -= m[k1+0,k2,k3] / δ[1]
    end
    for k1=1:n1, k2=1:n2-1, k3=1:n3
        d[k1,k2,k3,2] += m[k1,k2+1,k3] / δ[2]
        d[k1,k2,k3,2] -= m[k1,k2+0,k3] / δ[2]
    end
    for k1=1:n1, k2=1:n2, k3=1:n3-1
        d[k1,k2,k3,3] += m[k1,k2,k3+1] / δ[3]
        d[k1,k2,k3,3] -= m[k1,k2,k3+0] / δ[3]
    end
    d
end

function JopGradient_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,4}; δ, kwargs...) where {T}
    m .= 0
    n1,n2,n3 = size(m)
    for k1=1:n1-1, k2=1:n2, k3=1:n3
        m[k1+1,k2,k3] += d[k1,k2,k3,1] / δ[1]
        m[k1+0,k2,k3] -= d[k1,k2,k3,1] / δ[1]
    end
    for k1=1:n1, k2=1:n2-1, k3=1:n3
        m[k1,k2+1,k3] += d[k1,k2,k3,2] / δ[2]
        m[k1,k2+0,k3] -= d[k1,k2,k3,2] / δ[2]
    end
    for k1=1:n1, k2=1:n2, k3=1:n3-1
        m[k1,k2,k3+1] += d[k1,k2,k3,3] / δ[3]
        m[k1,k2,k3+0] -= d[k1,k2,k3,3] / δ[3]
    end
    m
end
