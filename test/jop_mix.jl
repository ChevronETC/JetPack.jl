using JetPack, Jets, Test

n1,n2 = 101,201

@testset "JopMix, linearity test, T=$T, mix=$mix" for T in (Float32,Float64), mix in ((4,3), (0,3), (3,0))
    A = JopMix(JetSpace(T,n1,n2), mix)
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopMix, dot product test, T=$T, mix=$mix" for T in (Float32,Float64), mix in ((4,3), (0,3), (3,0))
    A = JopMix(JetSpace(T,n1,n2), mix)
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
