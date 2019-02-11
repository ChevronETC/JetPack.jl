using JetPack, Jets, LinearAlgebra, Test

@testset "JopDifference, linearity, n=$n, T=$T" for n in ( (11,), (11,22), (11,22,33) ), T in (Float32,Float64)
    A = JopDifference(JetSpace(T,n))
    lhs, rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopDifference, dot-product, n=$n, T=$T" for n in ( (11,), (11,22), (11,22,33) ), T in (Float32,Float64)
    A = JopDifference(JetSpace(T,n))
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopDifference, correctness, n=$n" for n in ( (11,), (11,12), (11,12,13) )
    ndims = length(n)
    ops = [JopDifference(JetSpace(Float64,n), i) for i=1:ndims]
    α = [rand(Float64) for i=1:ndims]
    x = [mapreduce(i->α[i]*I.I[i], +, 1:length(n)) for I in CartesianIndices(n)]
    y = [A*x for A in ops]
    δ = [abs(sum(y[i] .- α[i]))/length(y[i]) for i=1:ndims]
    map(i->@test(δ[i]<0.1), 1:length(n))
end
