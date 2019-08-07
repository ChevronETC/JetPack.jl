"""
    op = JopReal(dom)

Extract the real part of a complex input array where `dom::JetSpace{<:Complex}`.
"""
JopReal(dom::JetSpace{Complex{T}}) where {T} = JopLn(dom = dom, rng = JetSpace(T,size(dom)), df! = JopReal_df!, df′! = JopReal_df′!)
export JopReal

JopReal_df!(d::AbstractArray{T}, m::AbstractArray{Complex{T}}; kwargs...) where {T} = d .= real.(m)
JopReal_df′!(m::AbstractArray{Complex{T}}, d::AbstractArray{T}; kwargs...) where {T} = m .= d .+ 0 * im
