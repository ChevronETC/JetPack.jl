using LinearAlgebra, Jets, JetPack, SpecialFunctions, Test

n1,n2 = 33,44

@testset "JopErf, correctness T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    x1 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    @test F*x1 ≈ erf.(x1)
end

@testset "JopErf, odd funcion T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    x1 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    @test F * (-1 .* x1) ≈ -(F * x1)
end

@testset "JopErf, conjugate T=$(T)" for T in (Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    x1 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    @test erf.(conj.(x1)) ≈ conj.(erf.(x1))
end

@testset "JopErf, linearity test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    m0 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    J  = jacobian!(F, m0)
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
end

@testset "JopErf, dot product test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    m0 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    J  = jacobian!(F, m0)
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

@testset "JopErf, linearization test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopErf(JetSpace(T,n1,n2))
    m0 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    δm = -1 .+ 2 .* rand(domain(F))
    μ = sqrt.([1.0,1.0/2.0,1.0/4.0,1.0/8.0,1.0/16.0,1.0/32.0,1.0/64.0,1.0/128.0,1.0/256.0,1.0/512.0])
    observed, expected = linearization_test(F, m0, μ = μ, δm = δm)
    @show observed
    @show expected
    δ = minimum(abs, observed .- expected)
    @test δ < 0.1
end
