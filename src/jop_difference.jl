"""
    A = JopDifference(R::JetSpace[, dim=1])

`A*m` is similar to `diff(a,dims=dim)` where `a` is a vector in the space `R`.
In other words, it is the one-sided difference of `a` along the dimension `dim`.
"""
function JopDifference(spc::JetAbstractSpace, dim::Integer=1)
    JopLn(df! = JopDifference_df!, df′! = JopDifference_df′!, dom = spc, rng = spc, s = (dim=dim,))
end
export JopDifference

function JopDifference_df!(d, m; dim, kwargs...)
    d .= 0
    n,Rpre,Rpost = JopDifference_indices(m, dim)
    JopDifference_df_dim!(d, m, n, Rpre, Rpost)
end

function JopDifference_df′!(m, d; dim, kwargs...)
    m .= 0
    n,Rpre,Rpost = JopDifference_indices(m, dim)
    JopDifference_df′_dim!(m, d, n, Rpre, Rpost)
end

function JopDifference_indices(m, dim)
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)
    n,Rpre, Rpost
end

function JopDifference_df_dim!(d, m, n, Rpre, Rpost)
    for Ipost in Rpost
        for i = 1:n-1
            for Ipre in Rpre
                d[Ipre, i, Ipost] += m[Ipre, i+1, Ipost]
                d[Ipre, i, Ipost] -= m[Ipre, i, Ipost]
            end
        end
    end
    d
end

function JopDifference_df′_dim!(m, d, n, Rpre, Rpost)
    for Ipost in Rpost
        for i = 1:n-1
            for Ipre in Rpre
                m[Ipre, i+1, Ipost] += d[Ipre, i, Ipost]
                m[Ipre, i, Ipost] -= d[Ipre, i, Ipost]
            end
        end
    end
    m
end
