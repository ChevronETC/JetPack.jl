"""
    A = JopLeakyIntegration(spc; α=0.99, dim=1)

`A` is Leaky integration (i.e. a running sum) along dimension `dim`
of `spc::JetSpace{T,N}`.  The result is scaled by `α`.
"""
function JopLeakyIntegration(spc::JetSpace{T}; α=0.99, dim=1) where {T}
    JopLn(dom = spc, rng = spc, df! = JopLeakyIntegration_df!, df′! = JopLeakyIntegration_df′!,
        s = (dim=dim, α=T(α)))
end

export JopLeakyIntegration

JopLeakyIntegration_df!(d::AbstractArray, m::AbstractArray; dim, α, kwargs...) = leaky_helper!(d, m, dim, α, false)
JopLeakyIntegration_df′!(m::AbstractArray, d::AbstractArray; dim, α, kwargs...) = leaky_helper!(m, d, dim, α, true)

function leaky_helper!(d, m, dim, α, isadj)
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)
    N = isadj ? (n:-1:1) : (1:n)
    @show N
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
