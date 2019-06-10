using JetPack, Jets, Test

@testset "JopLaplacian, T=$(T)" for T in (Float32, Float64)
    T=Float64
    sp = JetSpace(T,64,128)
    A = JopLaplacian(sp)

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs

    sp = JetSpace(T,16,32,64)
    A = JopLaplacian(sp)

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
