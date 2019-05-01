function JopTranslation(p1::Array{T,2}, p2::Array{T,2}) where {T}
    n = size(p1)
    @assert size(p1) == size(p2)
    JopLn(dom = JetSpace(T,n), rng = JetSpace(T,n), df! = JopTranslation_df!, df′! = JopTranslation_df′!, s = (p1=p1, p2=p2))
end

export JopTranslation

JopTranslation_df!(d, m; p1, p2, kwargs...) = JopTranslation_helper!(d, m, p1, p2, false)
JopTranslation_df′!(m, d; p1, p2, kwargs...) = JopTranslation_helper!(m, d, p1, p2, true)

function JopTranslation_helper!(d::AbstractArray{T,2}, m::AbstractArray{T,2}, p1, p2, adj) where {T}
    n1,n2 = size(d)

    d .= 0
    for i1 = 1:n1, i2 = 1:n2
        x1 = i1 + p1[i1,i2]
        x2 = i2 + p2[i1,i2]

        j1 = floor(Int, x1)
        j2 = floor(Int, x2)

        if j1 ∉ 1:n1-1 || j2 ∉ 1:n2-1
            continue
        end

        w1 = (j1 + 1 - x1) * (j2 + 1 - x2)
        w2 = (x1 - j1)     * (j2 + 1 - x2)
        w3 = (j1 + 1 - x1) * (x2 - j2)
        w4 = (x1 - j1)     * (x2 - j2)

        if adj
            d[j1,  j2]   += w1*m[i1,i2]
            d[j1+1,j2]   += w2*m[i1,i2]
            d[j1,  j2+1] += w3*m[i1,i2]
            d[j1+1,j2+1] += w4*m[i1,i2]
        else
            d[i1,i2] = w1*m[j1,j2] + w2*m[j1+1,j2] + w3*m[j1,j2+1] + w4*m[j1+1,j2+1]
        end
    end
    d
end
