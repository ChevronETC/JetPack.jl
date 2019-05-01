using JetPack, Jets, Test

@testset "JopPermutedims, 2D" begin
    A = JopPermutedims(JetSpace(Float64,3,2),[2;1])
    m = rand(domain(A))
    d = A*m
    @test d ≈ transpose(m)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopPermutedims, 3D" begin
    A = JopPermutedims(JetSpace(Float64,4,3,2),[3;1;2])
    m = rand(domain(A))
    d = A*m
    @test d ≈ permutedims(m,[3;1;2])
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
