"""
	op = JopImag(dom)

Extract the imaginary part of a complex input array where `dom::JetSpace{<:Complex}`.
"""

JopImag(dom::JetSpace{Complex{T}}) where {T} = JopLn(dom = dom, rng = JetSpace(T,size(dom)), df! = JopImag_df!, df′! = JopImag_df′!)
export JopImag

JopImag_df!(d::AbstractArray{T}, m::AbstractArray{Complex{T}}; kwargs...) where {T} = d .= imag.(m)
JopImag_df′!(m::AbstractArray{Complex{T}}, d::AbstractArray{T}; kwargs...) where {T} = m .= im .* d
