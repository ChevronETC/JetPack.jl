function JopProjection(u::AbstractArray{T,N}) where {T,N}
    u_normed = zeros(eltype(u), size(u))
    for idx in CartesianIndices(size(u)[1:end-1])
        u_normed[idx,:] .= u[idx,:] ./ norm(u[idx,:])
    end
    n = size(u)
    dom = JetSpace(T, n)
    rng = JetSpace(T, n[1:end-1])
    JopLn(dom = dom, rng = rng, df! = JopProjection_df!, df′! = JopProjection_df′!, s = (u=u_normed,))
end

export JopProjection

function JopProjection_df!(d::AbstractArray, m::AbstractArray; u, kwargs...)
    for idx in CartesianIndices(size(d))
        d[idx] = dot(u[idx,:], m[idx,:])
    end
    d
end

function JopProjection_df′!(m::AbstractArray, d::AbstractArray; u, kwargs...)
    for idx in CartesianIndices(size(d))
        for i = 1:size(m,ndims(m))
            m[idx,i] = u[idx,i]*d[idx]
        end
    end
    m
end
