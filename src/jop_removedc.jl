"""
    JopReoveDC(s::JetAbstractSpace{T,N})

F(x) = x - mean(x)
```
"""
function JopRemoveDC(spc::JetAbstractSpace{T,N}) where {T,N}
    JopLn(dom = spc, rng = spc, df! = JopRemoveDC_df!, df′! = JopRemoveDC_df′!)
end

export JopRemoveDC

function JopRemoveDC_df!(rngvec::AbstractArray{T,2}, domvec::AbstractArray{T,2}; kwargs...) where {T}
    rngvec .= domvec .- Statistics.mean(domvec)
end

function JopRemoveDC_df′!(domvec::AbstractArray{T,2}, rngvec::AbstractArray{T,2}; kwargs...) where {T}
    domvec .= rngvec .- Statistics.mean(rngvec)
end
