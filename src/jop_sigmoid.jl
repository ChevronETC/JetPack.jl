"""
    F = JopSigmoid(spc, c)

where `F` is the sigmoid operator 1/(1+e^(-c/x)) with domain and range given by `spc::JetSpace`, and scalar value c.
"""
function JopSigmoid(spc::JetSpace{T}, c = 1) where {T} 
    tmp = zeros(spc)
    JopNl(dom = spc, rng = spc, f! = JopSigmoid_f!, df! = JopSigmoid_df!, df′! = JopSigmoid_df′!, s = (c=T(c), tmp=tmp))
end
export JopSigmoid

JopSigmoid_f!(d::AbstractArray{T}, m::AbstractArray{T}; c, tmp) where {T} = d .= (T(1) .+ exp.(-m ./ c)).^T(-1)

function JopSigmoid_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ, c, tmp) where {T}
    tmp .= ( T(1) .+ exp.(-mₒ ./ c) ).^T(-1)
    δd .= δm .* (T(1) ./ c) .* (tmp .* (T(1) .- tmp))
end

function JopSigmoid_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ, c, tmp) where {T}
    tmp .= ( T(1) .+ exp.(-mₒ ./ c) ).^T(-1)
    δm .= δd .* conj((T(1) ./ c) .* (tmp .* (T(1) .- tmp)))
end
