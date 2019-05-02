using JetPack, Jets, Test

@testset "JopLMO" begin
    A = JopLMO(JetSpace(Float64,1024,256))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
