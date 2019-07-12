using JetPack, Jets, Test

n1,n2 = 101,201

@testset "JopRemoveDC, correctness test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64}) 
    A = JopRemoveDC(JetSpace(T,n1,n2))
    x = -1 .+ 2 .* rand(domain(A)) .+ rand(T)
    y1 = x .- (sum(x) / length(x))
    y2 = A * x
    @test y1 ≈ y2
end

@testset "JopRemoveDC, linearity test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64}) 
    A = JopRemoveDC(JetSpace(T,n1,n2))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopRemoveDC, dot product test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64})
    A = JopRemoveDC(JetSpace(T,n1,n2))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
