using LinearAlgebra, Jets, JetPack, Test, Flux

spc = JetSpace(Float64, 10, 10)

### set up a nonlinear operator according to the function `f`
f(x) = x ^ 2
F = JopAD(f; dom=spc, rng=spc)

@testset "JopAD: linearity, dot product, and linearization tests" begin
    
    x = rand(domain(F))
    J = jacobian(F, x)
    lhs, rhs = dot_product_test(J, rand(domain(J)), rand(range(J)))
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
    m0 = rand(domain(F))
    μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096,1/8192])
    δm = rand(domain(F))
    observed, expected = linearization_test(F, m0, μ = μ, δm = δm)
    δ = minimum(abs, observed - expected)
    @test δ < 1e-6
    
end