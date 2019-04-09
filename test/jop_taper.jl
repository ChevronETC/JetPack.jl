using JetPack, Jets, Test

@testset "taper" begin
    A = JopTaper(JetSpace(Float64,128,256), (2,), (.5,))
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test isapprox((lhs-rhs)/(lhs+rhs), 0.0, atol=1e-8)

    A = JopTaper(JetSpace(Float64,128,256), (2,), (.5,), mode=:fft)
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test isapprox((lhs-rhs)/(lhs+rhs), 0.0, atol=1e-8)

    A = JopTaper(JetSpace(Float64,128,256), (2,), (0.0,), (0.2,))
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test isapprox((lhs-rhs)/(lhs+rhs), 0.0, atol=1e-8)

    A = JopTaper(JetSpace(Float64,128,256), (2,), (0.0,), (0.2,), mode=:fft)
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test isapprox((lhs-rhs)/(lhs+rhs), 0.0, atol=1e-8)
end
