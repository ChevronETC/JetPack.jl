"""
    F = JopErf(spc)

where `F` is the error function operator with domain and range given by `spc::JetSpace`.
we use 'erf'(z) = 2/√π ∫_{0}^{z} \exp{-t^2}dt,
the derivative of which is 2/√π  \exp{-z^2}.
We expect the domain and range to be real.
"""
JopErf(spc::JetSpace) = JopNl(dom = spc, rng = spc, f! = JopErf_f!, df! = JopErf_df!)
export JopErf

# TODO: assert denominator in linearization is > 0

JopErf_f!(d::AbstractArray{T}, m::AbstractArray{T}) where {T<:Real} = d .= erf.(m)
JopErf_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ) where {T<:Real} = δd .= (2/sqrt( π )) .* exp.(-mₒ.^2) .* δm