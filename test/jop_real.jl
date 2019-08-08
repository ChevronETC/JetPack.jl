using Jets, JetPack, LinearAlgebra, Test

n1,n2 = 33,44

@testset "JopReal, dot product test, T=$(T)" for T in (Float64,Float32)
    A = JopReal(JetSpace(Complex{T},n1,n2))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopReal, linearity test, T=$(T)" for T in (Float64,Float32)
    A = JopReal(JetSpace(Complex{T},n1,n2))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopReal, correctness test, T=$(T)" for T in (Float64,Float32)
    A = JopReal(JetSpace(Complex{T},n1,n2))
    x = rand(domain(A))
    a = A * x
    b = real.(x)
    @test a ≈ b
end
