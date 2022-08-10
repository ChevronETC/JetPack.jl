"""
    F = JopAtan(spc, c)

where `F` is the shifted arc tangent operator F = arctan( x(t)/c ) + π/2 

x(t) = 0 corresponds to F = π/2
if |x(t)| ≤ c, then π/4 ≤ F ≤ 3π/4

    my change
"""
JopAtan(spc::JetSpace{T}, c = 1) where {T} = JopNl(dom = spc, rng = spc, f! = JopAtan_f!, df! = JopAtan_df!, df′! = JopAtan_df′!, s = (c=T(c),))
export JopAtan

function JopAtan_f!(d::AbstractArray, m::AbstractArray; c)
    d .= atan.(m ./ c) .+ (π/2)
end

function JopAtan_df!(δd::AbstractArray, δm::AbstractArray; mₒ, c) 
    δd .= (c ./ (c^2 .+ mₒ.^2)) .* δm
end

function JopAtan_df′!(δm::AbstractArray, δd::AbstractArray; mₒ, c) 
    δm .= conj(c ./ (c^2 .+ mₒ.^2)) .* δd
end
