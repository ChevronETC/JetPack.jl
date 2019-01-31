using JetPack, Jets, Test

@testset "JopBlend" begin
    A = JopBlend(Float64, 128, [5,10,20])
    lhs,rhs = dot_product_test(A,rand(domain(A)),rand(range(A)))
    @test lhs ≈ rhs
end
