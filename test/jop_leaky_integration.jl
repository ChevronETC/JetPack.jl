using JetPack, Jets, Test

@testset "JopLeakyIntegration, 2D, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64})
    nz, nx = 10,11
    A = JopLeakyIntegration(JetSpace(T,nz,nx), α=1.0)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs

    m = ones(domain(A))
    d = A*m
    @test real(d[1,1]) ≈ 1.0
    @test real(d[end,1]) ≈ nz
    @test real(d[end,end]) ≈ nz
end

@testset "JopLeakyIntegration, 3D, T=$T" for T in (Float32, Float64, Complex{Float32}, Complex{Float64})
    nz, ny, nx = 10,11,12
    A = JopLeakyIntegration(JetSpace(T,nz,ny,nx), α=1.0)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs

    m = ones(domain(A))
    d = A*m
    @test real(d[1,1,1]) ≈ 1.0
    @test real(d[end,1,1]) ≈ nz
    @test real(d[end,1,end]) ≈ nz
    @test real(d[end,end,1]) ≈ nz
    @test real(d[end,end,end]) ≈ nz
end
