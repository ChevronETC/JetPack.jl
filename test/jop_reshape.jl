using JetPack, Jets, Test

n1,n2 = 101,201

@testset "JopReshape, linearity test, T=$T" for T in (Float64,Float32)
    A = JopReshape(JetSpace(T,n1,n2), JetSpace(T,n1,n2,1))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopReshape, dot product test, T=$T" for T in (Float64,Float32)
    A = JopReshape(JetSpace(T,n1,n2,1), JetSpace(T,n1,n2))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs rtol=1e-6

    A = JopReshape(JetSpace(T,n1,n2), JetSpace(T,n1,n2,1))
    @test lhs ≈ rhs rtol=1e-6
end

@testset "JopReshape, correctness" begin
    A = JopReshape(JetSpace(Float64,n1*n2), JetSpace(Float64,n1,n2))
    m = rand(domain(A))
    d = A * m
    @test isa(d, Array{Float64,2})
    @test d[:] ≈ m[:]
    a = A' * d
    @test isa(a, Vector{Float64})
    @test a[:] ≈ m[:]
end
