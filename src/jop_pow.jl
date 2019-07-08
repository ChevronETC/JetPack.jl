"""
    F = JopPow(spc, c, a)

where `F` is the power operator (x/c)^a with domain and range given by `spc::JetSpace`, and scalar values c, a.
"""
JopPow(spc::JetSpace{T}, c = 1, a = 2) where {T} = JopNl(dom = spc, rng = spc, f! = JopPow_f!, df! = JopPow_df!, df′! = JopPow_df′!, s = (c=T(c), a=T(a)))
export JopPow

JopPow_f!(d::AbstractArray, m::AbstractArray; c, a) = d .= (m ./ c).^a
JopPow_df!(δd::AbstractArray, δm::AbstractArray; mₒ, c, a) = δd .= ((a/c) .* (mₒ ./ c).^(a-1)) .* δm
JopPow_df′!(δm::AbstractArray, δd::AbstractArray; mₒ, c, a)= δm .= conj((a/c) .* (mₒ ./ c).^(a-1)) .* δd
