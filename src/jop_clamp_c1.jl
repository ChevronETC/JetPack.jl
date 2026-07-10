"""
    F = JopClampC1(spc(T), lo::T, hi::T; transition=0.05)

Pointwise C1 clamp operator with domain and range given by `spc::JetSpace`.
The operator acts like the identity in the interior of `[lo, hi]`, saturates to the
hard bounds outside that interval, and uses a cubic transition of width
`δ = transition * (hi - lo)` near each bound to keep the map continuously
differentiable.

More precisely, the lower transition on `[lo, lo + δ]` is

`lo + δ * (-s^3 + 2s^2)`, with `s = (x - lo) / δ`,

and the upper transition on `[hi - δ, hi]` is defined symmetrically.

This operator is useful when hard clamping would create derivative jumps that make
linearization tests or gradient-based methods unnecessarily brittle.
"""
function JopClampC1(spc::JetSpace{T}, lo::Real, hi::Real; transition::Real=0.05) where {T<:AbstractFloat}
    @assert hi > lo "JopClampC1 -- require hi > lo. Got lo=$(lo), hi=$(hi)."
    _transition = T(transition)
    @assert zero(T) <= _transition <= one(T)/2 "JopClampC1 -- require transition in [0, 0.5]. Got transition=$transition."
    JopNl(dom = spc, rng = spc, f! = JopClampC1_f!, df! = JopClampC1_df!, df′! = JopClampC1_df′!, s = (lo=T(lo), hi=T(hi), transition=_transition))
end
export JopClampC1

@inline function _clamp_c1_forward(x::T, lo::T, hi::T, δ::T) where {T<:AbstractFloat}
    if x <= lo
        return lo
    elseif x < lo + δ
        s = (x - lo) / δ
        return lo + δ * (-s*s*s + 2*s*s)
    elseif x <= hi - δ
        return x
    elseif x < hi
        s = (hi - x) / δ
        return hi - δ * (-s*s*s + 2*s*s)
    else
        return hi
    end
end

@inline function _clamp_c1_derivative(x::T, lo::T, hi::T, δ::T) where {T<:AbstractFloat}
    if x <= lo
        return zero(T)
    elseif x < lo + δ
        s = (x - lo) / δ
        return -3*s*s + 4*s
    elseif x <= hi - δ
        return one(T)
    elseif x < hi
        s = (hi - x) / δ
        return -3*s*s + 4*s
    else
        return zero(T)
    end
end

function JopClampC1_f!(d::AbstractArray{T}, m::AbstractArray{T}; lo, hi, transition) where {T<:AbstractFloat}
    δ = transition * (hi - lo)
    if δ <= eps(T)
        d .= clamp.(m, lo, hi)
        return d
    end

    @inbounds @simd for i in eachindex(d, m)
        d[i] = _clamp_c1_forward(m[i], lo, hi, δ)
    end
    d
end

function JopClampC1_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ, lo, hi, transition) where {T}
    δ = transition * (hi - lo)
    if δ <= eps(T)
        @inbounds @simd for i in eachindex(δd, δm, mₒ)
            inside = (mₒ[i] >= lo) & (mₒ[i] <= hi)
            δd[i] = inside ? δm[i] : zero(T)
        end
        return δd
    end

    @inbounds @simd for i in eachindex(δd, δm, mₒ)
        δd[i] = _clamp_c1_derivative(mₒ[i], lo, hi, δ) * δm[i]
    end
    δd
end

function JopClampC1_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ, lo, hi, transition) where {T}
    δ = transition * (hi - lo)
    if δ <= eps(T)
        @inbounds @simd for i in eachindex(δm, δd, mₒ)
            inside = (mₒ[i] >= lo) & (mₒ[i] <= hi)
            δm[i] = inside ? conj(one(T)) * δd[i] : zero(T)
        end
        return δm
    end

    @inbounds @simd for i in eachindex(δm, δd, mₒ)
        δm[i] = conj(_clamp_c1_derivative(mₒ[i], lo, hi, δ)) * δd[i]
    end
    δm
end
