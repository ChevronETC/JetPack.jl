using JetPack, Jets, Test, LinearAlgebra

n1,n2,n3 = 10,11,12

@testset "JopLeakyIntegration, 1D, correctness, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    n1 = 6
    A = JopLeakyIntegration(JetSpace(T,n1), α=α)
    x = ones(domain(A))
    y1 = A * x
    y2 = zeros(domain(A))
    y2[1] = T(1   + 0   + 0   + 0   + 0   +   0)
    y2[2] = T(α^1 + 1   + 0   + 0   + 0   +   0)
    y2[3] = T(α^2 + α^1 + 1   + 0   + 0   +   0)
    y2[4] = T(α^3 + α^2 + α^1 + 1   + 0   +   0)
    y2[5] = T(α^4 + α^3 + α^2 + α^1 + 1   +   0)
    y2[6] = T(α^5 + α^4 + α^3 + α^2 + α^1 +   1)
    @test norm(y1 .- y2) < 10 * eps(real(T))
end

@testset "JopLeakyIntegration, 1D, adjoint, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1), α=α)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 2D, adjoint, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1,n2), α=α)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 3D, adjoint, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1,n2,n3), α=α)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 1D, linearity, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1), α=α)
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(A')
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 2D, linearity, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1,n2), α=α)
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(A')
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 2D, linearity, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1,n2), α=α)
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(A')
    @test lhs ≈ rhs
end

@testset "JopLeakyIntegration, 3D, linearity, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64}), α in (1.0, real(rand(T)))
    A = JopLeakyIntegration(JetSpace(T,n1,n2,n3), α=α)
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(A')
    @test lhs ≈ rhs
end
