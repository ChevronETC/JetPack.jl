"""
```op = JopImag(dom,rng)```

This operator takes the imaginary part of a required complex domain input and returns a real range
	dom::JetSpace{Complex{T},N}` is the domain of the operator.
	rng::JetSpace{T,N}` is the range of the operator.

im A = 1/(2i) (A - A*)
"""

function JopImag(dom::JetSpace{Complex{T},N}) where {T,N}
    JopLn(dom = dom, rng = JetSpace(T,size(dom)), df! = JopImag_df!, df′! = JopImag_df′!)
end
export JopImag

function JopImag_df!(rngvec::AbstractArray{T,N}, domvec::AbstractArray{Complex{T},N}; kwargs...) where {T,N}
	rngvec .= imag.(domvec)
end

function JopImag_df′!(domvec::AbstractArray{Complex{T},N}, rngvec::AbstractArray{T,N}; kwargs...) where {T,N}
	domvec .= im .* rngvec
end
