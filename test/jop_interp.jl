using JetPack, Jets, Test

n1dom,n2dom = 10,11
n1rng,n2rng = 29,31

@testset "JopInterp, linearity test" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom,n2dom), JetSpace(T,n1rng,n2rng))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopInterp, dot product test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom,n2dom), JetSpace(T,n1rng,n2rng))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopInter, test bugfix for specific dom,rng, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,280,222), JetSpace(T,281,2221))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
