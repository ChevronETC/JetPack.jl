using JetPack, Jets, Test

@testset "JopProjection, 2D, first component" begin
    u = zeros(256,2)
    u[:,1] .= 1.0
    u[:,2] .= 0.0

    P = JopProjection(u)
    m = zeros(256,2)
    m[:,1] .= 1.0
    m[:,2] .= 2.0

    d = P*m
    a = P'*d

    d_expected = 1.0*ones(256)
    a_expected = zeros(256,2)
    a_expected[:,1] .= 1.0

    @test d ≈ 1.0*ones(256)
    @test a ≈ a_expected

    lhs,rhs = dot_product_test(P, rand(domain(P)), rand(range(P)))
    @test lhs ≈ rhs
end

@testset "JopProjection, 2D, second component" begin
    u = zeros(256,2)
    u[:,1] .= 0.0
    u[:,2] .= 1.0

    P = JopProjection(u)
    m = zeros(256,2)
    m[:,1] .= 1.0
    m[:,2] .= 2.0

    d = P*m
    a = P'*d

    d_expected = 2.0*ones(256)
    a_expected = zeros(256,2)
    a_expected[:,2] .= 2.0

    @test d ≈ 2.0*ones(256)
    @test a ≈ a_expected

    lhs,rhs = dot_product_test(P, rand(domain(P)), rand(range(P)))
    @test lhs ≈ rhs
end
