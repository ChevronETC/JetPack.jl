using JetPack, Jets, Test

@testset "highpass" for T in (Float32, Float64)
    A = JopHighpass(JetSpace(T,64,128))

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs rtol=1e-5
end
