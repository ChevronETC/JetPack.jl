using LinearAlgebra, Jets, JetPack, Test

@inline function _clamp_c1_forward_test(x::T, lo::T, hi::T, δ::T) where {T<:AbstractFloat}
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

n1, n2 = 17, 9

function _mixed_state(T, dims, lo, hi; transition=T(0))
    m0 = fill((lo + hi) / 2, dims...)
    if transition > 0
        δ = transition * (hi - lo)
        @views m0[1:3, :] .= lo - T(0.30) * δ
        @views m0[4:6, :] .= lo + T(0.30) * δ
        @views m0[7:9, :] .= (lo + hi) / 2
        @views m0[10:12, :] .= hi - T(0.30) * δ
        @views m0[13:15, :] .= hi + T(0.30) * δ
    else
        gap = T(0.4)
        @views m0[1:5, :] .= lo - gap
        @views m0[6:11, :] .= (lo + hi) / 2
        @views m0[12:17, :] .= hi + gap
    end
    m0
end

function _oscillatory_step(T, dims; scale=T(0.08))
    δm = Array{T}(undef, dims...)
    @inbounds for i in eachindex(δm)
        δm[i] = scale * sin(T(i))
    end
    δm
end

@testset "JopClampC1, correctness T=$(T), transition=$(transition)" for T in (Float32, Float64), transition in (0.0, 0.05)
    lo = T(1.5)
    hi = T(5.5)
    F = JopClampC1(JetSpace(T, 7), lo, hi; transition=transition)
    x = T[0.0, lo, lo + T(0.1), (lo + hi) / 2, hi - T(0.1), hi, hi + T(1.0)]
    y = F * x

    if transition == 0
        expected = clamp.(x, lo, hi)
    else
        δ = T(transition) * (hi - lo)
        expected = similar(x)
        @inbounds for i in eachindex(x)
            expected[i] = _clamp_c1_forward_test(x[i], lo, hi, δ)
        end
    end

    @test y ≈ expected
end

@testset "JopClampC1, linearity and dot product T=$(T), transition=$(transition)" for T in (Float32, Float64), transition in (0.0, 0.05)
    lo = T(1.5)
    hi = T(5.5)
    transition = T(transition)
    F = JopClampC1(JetSpace(T, n1, n2), lo, hi; transition=transition)

    m0 = _mixed_state(T, size(domain(F)), lo, hi; transition=transition)

    J = jacobian(F, m0)
    lhs, rhs = linearity_test(J)
    @test lhs ≈ rhs

    J = jacobian!(F, m0)
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

@testset "JopClampC1, linearization rates T=$(T), transition=$(transition)" for T in (Float32, Float64), transition in (0.0, 0.05)
    lo = T(1.5)
    hi = T(5.5)
    transition = T(transition)
    F = JopClampC1(JetSpace(T, n1, n2), lo, hi; transition=transition)

    if transition == 0
        m0 = _mixed_state(T, size(domain(F)), lo, hi; transition=transition)
        δm = _oscillatory_step(T, size(m0); scale=T(0.08))
    else
        m0 = _mixed_state(T, size(domain(F)), lo, hi; transition=transition)
        δm = _oscillatory_step(T, size(m0))
    end

    J = jacobian!(F, m0)
    δd = J * δm

    μs = T[1e-1, 5e-2, 2.5e-2, 1.25e-2]
    e0 = T[norm(F * (m0 .+ μ .* δm) - F * m0) for μ in μs]
    e1 = T[norm(F * (m0 .+ μ .* δm) - F * m0 - μ .* δd) for μ in μs]

    rate0 = log(e0[1] / e0[end]) / log(μs[1] / μs[end])
    rate1 = log(e1[1] / e1[end]) / log(μs[1] / μs[end])

    if transition == 0
        @test abs(rate0 - 1) < 0.05
        @test maximum(e1) < (T === Float32 ? 1f-5 : 1e-12)
    else
        @test abs(rate0 - 1) < 0.05
        @test abs(rate1 - 2) < 0.05
    end
end

@testset "JopClampC1, hard clamp loses smooth rates at boundary T=$(T)" for T in (Float32, Float64)
    lo = T(1.5)
    hi = T(5.5)
    F = JopClampC1(JetSpace(T, n1, n2), lo, hi; transition=0)

    m0 = fill((lo + hi) / 2, size(domain(F))...)
    @views m0[1:6, :] .= lo .- T(1e-3)
    @views m0[7:12, :] .= (lo + hi) / 2
    @views m0[13:17, :] .= hi .+ T(1e-3)

    δm = zeros(T, size(m0)...)
    @views δm[1:6, :] .= T(0.2)
    @views δm[7:12, :] .= _oscillatory_step(T, size(m0[7:12, :]); scale=T(0.05))
    @views δm[13:17, :] .= -T(0.2)

    J = jacobian!(F, m0)
    δd = J * δm

    μs = T[1e-1, 5e-2, 2.5e-2, 1.25e-2]
    e0 = T[norm(F * (m0 .+ μ .* δm) - F * m0) for μ in μs]
    e1 = T[norm(F * (m0 .+ μ .* δm) - F * m0 - μ .* δd) for μ in μs]

    rate0 = log(e0[1] / e0[end]) / log(μs[1] / μs[end])
    rate1 = log(e1[1] / e1[end]) / log(μs[1] / μs[end])

    @test_broken abs(rate0 - 1) < 0.05
    @test_broken abs(rate1 - 2) < 0.05
end
