"""
    A = JopLeakyIntegration(spc; α=0.99, dim=1)

`A` is Leaky integration (i.e. a running sum) along dimension `dim`
of `spc::JetSpace{T,N}`.  The result is scaled by `α`.
"""

#=
From stanford: http://sep.stanford.edu/sep/prof/geelib/leakint.lop
    tt = 0.0

    if (adj)
        do i =  n, 1, -1 {
            tt = rho*tt + yy(i)
            xx(i) += tt
        }
    else
        do i =  1, n {
            tt = rho*tt + xx(i)
            yy(i) += tt
        }
    }
=#

function JopLeakyIntegration(spc::JetSpace{T}; α=0.99, dim=1) where {T}
    JopLn(dom = spc, rng = spc, df! = JopLeakyIntegration_df!, df′! = JopLeakyIntegration_df′!,
        s = (dim=dim, α=T(α)))
end

export JopLeakyIntegration

function JopLeakyIntegration_df!(d::AbstractArray{T}, m::AbstractArray{T}; dim, α, kwargs...) where {T}
    n1 = size(d,1)
    n2 = div(length(d),n1)
    _d = reshape(d,n1,n2)
    _m = reshape(m,n1,n2)
    for k2 = 1:n2
        tmp = T(0)
        for k1 = 1:n1
            tmp = α * tmp + _m[k1,k2] 
            _d[k1,k2] = tmp
        end
    end
    _d
end

function JopLeakyIntegration_df′!(m::AbstractArray{T}, d::AbstractArray{T}; dim, α, kwargs...) where {T}
    n1 = size(d,1)
    n2 = div(length(d),n1)
    _d = reshape(d,n1,n2)
    _m = reshape(m,n1,n2)
    for k2 = 1:n2
        tmp = T(0)
        for k1 = n1:-1:1
            tmp = α * tmp + _d[k1,k2] 
            _m[k1,k2] = tmp
        end
    end
    _m
end

#=
I think this is broken ... wont the first sample have value α?

JopLeakyIntegration_df!(d::AbstractArray, m::AbstractArray; dim, α, kwargs...) = leaky_helper!(d, m, dim, α, false)
JopLeakyIntegration_df′!(m::AbstractArray, d::AbstractArray; dim, α, kwargs...) = leaky_helper!(m, d, dim, α, true)

function leaky_helper!(d, m, dim, α, isadj)
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)
    N = isadj ? (n:-1:1) : (1:n)
    leaky_helper_dim!(d, m, α, N, Rpre, Rpost)
end

function leaky_helper_dim!(d::AbstractArray{T}, m, α, N, Rpre, Rpost) where {T}
    for Ipost in Rpost
        σ = zero(T)
        for i in N
            for Ipre in Rpre
                σ += m[Ipre, i, Ipost]
                d[Ipre, i, Ipost] = α*σ
            end
        end
    end
    d
end
=#
