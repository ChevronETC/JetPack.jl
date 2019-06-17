"""
    F = JopExp(spc, c)

where `F` is the exponential operator e^(c x) with domain and range given by `spc::JetSpace`, and scalar value c.
"""
JopExp(spc::JetSpace{T}, c = 1) where {T} = JopNl(dom = spc, rng = spc, f! = JopExp_f!, df! = JopExp_df!, df′! = JopExp_df′!, s = (c=T(c),))
export JopExp

JopExp_f!(d::AbstractArray, m::AbstractArray; c) = d .= exp.(c .* m)
JopExp_df!(δd::AbstractArray, δm::AbstractArray; mₒ, c) = δd .= c .* δm .* exp.(c .* mₒ)
JopExp_df′!(δm::AbstractArray, δd::AbstractArray; mₒ, c) = δm .= conj(c) .* δd .* conj.(exp.(c .* mₒ))
