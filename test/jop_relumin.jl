using LinearAlgebra, Jets, JetPack, Test

# note no need for a linearity test as function is "linear" for x > 0 and constant for x <= 0 
# you do need special handling: strictly positive or negative domain and range vectors for tests

n1,n2 = 33,44

@testset "JopReluMin, correctness T=$(T)" for T in (Float64,Float32)
    F = JopReluMin(JetSpace(T,n1,n2))
    x1 = rand(domain(F))
    @test F*x1 ≈ min.(x1,0)
end

@testset "JopReluMin, linearity test, T=$(T)" for T in (Float64,Float32)
    F = JopReluMin(JetSpace(T,n1,n2))
    m1 = rand(domain(F)) .+ T(0.1)
    m2 = rand(domain(F)) .+ T(0.1)
    lhs,rhs = linearity_test(F, m1, m2)
    @test lhs ≈ rhs

    m1 = - rand(domain(F)) .- T(0.1)
    m2 = - rand(domain(F)) .- T(0.1)
    lhs,rhs = linearity_test(F, m1, m2)
    @test lhs ≈ rhs
end

@testset "JopReluMin, dot product test, T=$(T)" for T in (Float64,Float32)
    F = JopReluMin(JetSpace(T,n1,n2))
    δm = rand(domain(F)) .+ T(0.1)
    δd = rand(domain(F)) .+ T(0.1)
    lhs, rhs = dot_product_test(F, δd, δm)
    @test lhs ≈ rhs

    F = JopReluMin(JetSpace(T,n1,n2))
    δm = - rand(domain(F)) .- T(0.1)
    δd = - rand(domain(F)) .- T(0.1)
    lhs, rhs = dot_product_test(F, δd, δm)
    @test lhs ≈ rhs
end

