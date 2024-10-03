using JetPack, Jets, Test

@testset "JopSpreadStack, functions, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64})
    A = JopSpreadStack(JetSpace(T,11), JetSpace(T,11,5))
    x = rand(domain(A))
    y = A * x
    for k ∈ 1:size(y,2)
        @test x[:] ≈ y[:,k]
    end
    z = A' * y
    @test z[:] ≈ size(y,2) * x[:]
end

@testset "JopSpreadStack, linearity test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64})
    A = JopSpreadStack(JetSpace(T,11), JetSpace(T,11,5))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopSpreadStack, dot product test, T=$T" for T in  (Float32,Float64,Complex{Float32},Complex{Float64})
    A = JopSpreadStack(JetSpace(T,11), JetSpace(T,11,5))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

nothing
