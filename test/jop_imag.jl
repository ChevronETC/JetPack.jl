using Jets, JetPack, LinearAlgebra, Test

n1,n2 = 33,44

@testset "JopImag, dot product test, T=$(T)" for T in (Float64,Float32)
    A = JopImag(JetSpace(Complex{T},n1,n2))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopImag, linearity test, T=$(T)" for T in (Float64,Float32)
    A = JopImag(JetSpace(Complex{T},n1,n2))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopImag, correctness test, T=$(T)" for T in (Float64,Float32)
    A = JopImag(JetSpace(Complex{T},n1,n2))
    x = rand(domain(A))
    a = A * x
    b = 1/(2*im) .* (x .- conj.(x))
    @test a ≈ b
end
