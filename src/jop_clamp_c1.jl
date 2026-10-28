"""
    F = JopClampC1(spc(T), lo::T, hi::T; transition=0.01)

Pointwise C1 clamp operator with domain and range given by `spc::JetSpace`.
The operator is the **exact identity on `[lo + δ, hi - δ]`**, where
`δ = transition * (hi - lo)`, so interior values are never modified. Outside that
band it bends smoothly toward the bounds using exponential transitions that are
**convex near `lo`** and **concave near `hi`**, approaching `lo`/`hi` asymptotically.

The lower transition on `x <= lo + δ` is

`f(x) = lo + δ * exp((x - lo - δ) / δ)`,

and the upper transition on `x >= hi - δ` is

`f(x) = hi - δ * exp((hi - δ - x) / δ)`.

Each transition is tangent to the identity at the band edge, so `f(x) >= x` just
below `lo + δ` and `f(x) <= x` just above `hi - δ`: values are nudged *up* toward the
lower bound and *down* toward the upper bound, never dragged across the identity. The
lower piece is convex (`f'' > 0`) and the upper piece concave (`f'' < 0`). The
derivative is `exp(...)` in the transitions (slope 1 at the band edge, decaying to 0 at
the asymptote) and `1` in the interior, so the map is C1. The bounds are approached
asymptotically and never exceeded. Passing `transition = 0` recovers a hard clamp.

This operator is useful when hard clamping would create derivative jumps that make
linearization tests or gradient-based methods unnecessarily brittle.
"""
function JopClampC1(spc::JetSpace{T}, lo::Real, hi::Real; transition::Real=0.01) where {T<:AbstractFloat}
    @assert hi > lo "JopClampC1 -- require hi > lo. Got lo=$(lo), hi=$(hi)."
    _transition = T(transition)
    @assert zero(T) <= _transition <= one(T)/2 "JopClampC1 -- require transition in [0, 0.5]. Got transition=$transition."
    JopNl(dom = spc, rng = spc, f! = JopClampC1_f!, df! = JopClampC1_df!, df′! = JopClampC1_df′!, s = (lo=T(lo), hi=T(hi), transition=_transition))
end
export JopClampC1

@inline function _clamp_c1_forward(x::T, lo::T, hi::T, δ::T) where {T<:AbstractFloat}
    if x <= lo + δ
        return lo + δ * exp((x - lo - δ) / δ)
    elseif x <= hi - δ
        return x
    else
        return hi - δ * exp((hi - δ - x) / δ)
    end
end

@inline function _clamp_c1_derivative(x::T, lo::T, hi::T, δ::T) where {T<:AbstractFloat}
    if x <= lo + δ
        return exp((x - lo - δ) / δ)
    elseif x <= hi - δ
        return one(T)
    else
        return exp((hi - δ - x) / δ)
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
