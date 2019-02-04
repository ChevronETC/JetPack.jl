using JetPack, Jets, Test

@testset "JopDiagonal, non constant diagonal, T=$T" for T in (Float64, Complex{Float64})
    diagonal = rand(T,512)
    A = JopDiagonal(diagonal)
    m = rand(domain(A))
    @test A*m ≈ diagonal .* m
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopDiagonal, non constant real diagonal, complex spaces" begin
    diagonal = rand(512)
    A = JopDiagonal(JetSpace(Complex{Float64}, 512), diagonal)
    m = rand(domain(A))
    @test eltype(m) <: Complex
    @test A*m ≈ diagonal .* m
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopDiagonal, constant diagonal" begin
    diagonal = rand(Float64)
    A = JopDiagonal(JetSpace(Float64,512), diagonal)
    m = rand(domain(A))
    @test A*m ≈ diagonal .* m
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
