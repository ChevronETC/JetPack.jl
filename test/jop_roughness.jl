using JetPack, Jets, Test

@testset "roughness, dot product test N=$N, dim=$dim" for (N,dim) in (
        ((256,512),1), ((16,32,64,4),1), ((16,32,64,4),2), ((16,32,64,4),3))
    A = JopRoughness(JetSpace(Float64,N),dim)

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "roughness, linearity test N=$N, dim=$dim" for (N,dim) in (
        ((256,512),1), ((16,32,64,4),1), ((16,32,64,4),2), ((16,32,64,4),3))
    A = JopRoughness(JetSpace(Float64,N),dim)

    lhs, rhs = linearity_test(A)
    @test lhs ≈ rhs
end
