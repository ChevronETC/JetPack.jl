"""
    F = JopTanh(spc)

where `F` is the hyperbolic tangent operator with domain and range given by `spc::JetSpace`.
we use 'tanh'(c*x) = (\exp{c*x}-\exp{-c*x})/(\exp{c*x}+\exp{-c*x}),
the derivative of which is 1-tanh(c*x)^2.
We expect the domain and range to be real.
"""
JopTanh(spc::JetSpace{T}, c = 1) where {T} = JopNl(dom = spc, rng = spc, f! = JopTanh_f!, df! = JopTanh_df!, df′! = JopTanh_df′!, s = (c=T(c),))
export JopTanh

JopTanh_f!(d::AbstractArray, m::AbstractArray; c) = d .= tanh.(c .* m)

JopTanh_df!(δd::AbstractArray, δm::AbstractArray; mₒ, c) = δd .= c .* (1 .- tanh.(c .* mₒ).^2) .* δm

JopTanh_df′!(δd::AbstractArray, δm::AbstractArray; mₒ, c) = δd .= conj.(c .* (1 .- tanh.(c .* mₒ).^2)) .* δm

