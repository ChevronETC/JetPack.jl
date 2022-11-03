"""
    A = JopBlend(T, nsamples, shottimes[, nsamples_blended])

where `A` is a shot mixing operation, so that if `d=Am`, then `m::Array{T,2}` is a
common receiver gather with each trace corresponding to a different shot.  The
excitation time of each are the entries of `shottimes::Array{Int,1}` in units of
time samples.  The size of `m` is `(nsamples,length(shottimes))`.  The size of
`d` is either `nsamples_blended` (if provided) or `maximum(shot_times)+nsamples-1`.

# Examples

## blend two receiver gather traces:
```julia
A = JopBlend(Float64, 128, [1,64])
m = rand(domain(A))
d = A*m # receiver gather with blended shots
```
"""
function JopBlend(T, nsamples, shottimes, nsamples_blended=0)
    if nsamples_blended == 0
        nsamples_blended = maximum(shottimes) + nsamples - 1
    end

    dom = JetSpace(T, nsamples, length(shottimes))
    rng = JetSpace(T, nsamples_blended)

    JopLn(df! = JopBlend_df!, df′! = JopBlend_df′!, dom = dom, rng = rng, s = (st=shottimes,))
end
export JopBlend

function JopBlend_df!(d::AbstractArray{T,1}, m::AbstractArray{T,2}; st, kwargs...) where {T}
    d .= 0
    nsamples, nshots = size(m)
    nsamples_blended = length(d)
    for ishot = 1:nshots, isample = 1:nsamples
        j = st[ishot] + isample - 1
        if j ∈ 1:nsamples_blended
            d[j] += m[isample,ishot]
        end
    end
    d
end

function JopBlend_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,1}; st, kwargs...) where {T}
    m .= 0
    nsamples, nshots = size(m)
    nsamples_blended = length(d)
    for ishot = 1:nshots, isample = 1:nsamples
        j = st[ishot] + isample - 1
        if j ∈ 1:nsamples_blended
            m[isample,ishot] = d[j]
        end
    end
    m
end
