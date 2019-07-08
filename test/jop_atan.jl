using LinearAlgebra, Jets, JetPack, Test, Printf

n1,n2 = 33,44

@testset "JopAtan, correctness T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopAtan(JetSpace(T,n1,n2))
    x1 = rand(domain(F))
    @test F*x1 ≈ atan.(x1) .+ π/2
end

@testset "JopAtan, linearity test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopAtan(JetSpace(T,n1,n2))
    J = jacobian(F, rand(domain(F)))
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
end

@testset "JopAtan, dot product test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopAtan(JetSpace(T,n1,n2))
    J  = jacobian!(F, rand(domain(F)))
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

@testset "JopAtan, linearization test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopAtan(JetSpace(T,n1,n2))
    m0 = 0.5 .* (-1 .+ rand(domain(F)))
    μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096,1/8192])
    δm = rand(domain(F))
    observed, expected = linearization_test(F, m0, μ = μ, δm = δm)
    δ = minimum(abs, observed - expected)
    @test δ < 0.1
end
