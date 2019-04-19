using JetPack, Jets, Test

@testset "restriction, real to real, 1D" begin
    A = JopRestriction(JetSpace(Float64,512), collect(1:4:512))
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
    m = rand(domain(A))
    @test A*m ≈ m[1:4:512]
    d = m[1:4:512]
    m = zeros(domain(A))
    m[1:4:512] .= d
    @test A'*d ≈ m
end

@testset "restriction, real to real, 2D" begin
    indices = [ (i,j) for i=[1,4], j=[1,4] ]
    A = JopRestriction(JetSpace(Float64,4,4), vec(indices))
    m = [11.0 12.0 13.0 14.0 ; 21.0 22.0 23.0 24.0 ; 31.0 32.0 33.0 34.0 ; 41.0 42.0 43.0 44.0]
    @test A*m ≈ [11.0, 41.0, 14, 44.0]
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
