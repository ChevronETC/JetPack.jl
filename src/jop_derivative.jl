"""
    A = JopDerivative(R::JetSpace[; dim=1, accuracy=4, delta=1.0)

`A*m` is the centered finite different approximation to the derivative of `m`
along dimension `dim`.  The `accuracy` of the approximation is either `4` for
a 4th order accurate estimate of the derivative, or `8` for an 8th order accurate
derivative.
"""
function JopDerivative(sp::JetSpace{T}; dim=1, accuracy=4, delta=1.0) where {T}
    local a
    if accuracy == 4
        a = T[1/12, -2/3, 0, 2/3, -1/12]
    elseif accuracy == 8
        a = T[1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]
    else
        error("order of accuracy must be either 4 or 8, got accuracy=$(accuracy)")
    end
    a ./= delta
    JopLn(df! = JopDerivative_df!, df′! = JopDerivative_df′!, dom = sp, rng = sp, s = (a=a,dim=dim))
end
export JopDerivative

JopDerivative_df!(d::AbstractArray, m::AbstractArray; a, dim, kwargs...) = derivative_helper(d, m, a, dim)
JopDerivative_df′!(m::AbstractArray, d::AbstractArray; a, dim, kwargs...) = derivative_helper(m, d, -a, dim)

function derivative_helper_dim(d, m, a, n, Rpre, Rpost)
    d .= 0
    M = div(length(a)-1,2)
    for Ipost in Rpost
        for i = (M+1):(n-M)
            for j = 1:length(a)
                k = i-M+j-1
                if M < k <= n-M
                    for Ipre in Rpre
                        d[Ipre, i, Ipost] += a[j] * m[Ipre, k, Ipost]
                    end
                end
            end
        end
    end
    d
end
function derivative_helper(d::AbstractArray, m::AbstractArray, a::Vector, dim::Integer)
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)
    derivative_helper_dim(d, m, a, n, Rpre, Rpost)
end
