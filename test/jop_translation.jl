using JetPack, Jets, Test

@testset "JopTranslation, dot product test" begin
    n1,n2 = 100,120
    A = JopTranslation(rand(n1,n2), rand(n1,n2))
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopTranslation, Linearity" begin
    n1,n2 = 100,120
    A = JopTranslation(randn(n1,n2), randn(n1,n2))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopTranslation, Correctness" begin
    n1,n2 = 100,120

    m = [1:n1;] * ones(n2)'
    A = JopTranslation(ones(n1,n2), zeros(n1,n2))
    d = A*m

    @test m[2:end-2,2:end-1] ≈ d[1:end-3,2:end-1]
end

@testset "JopTranslation, Correctness half width in t" begin
    n1,n2 = 100,120

    m = [1:n1;] * ones(n2)'
    A = JopTranslation(0.5*ones(n1,n2), 0.0*ones(n1,n2))
    d = A*m

    @test m[2:end-1,2:end-1] .+ 0.5 ≈ d[2:end-1,2:end-1]
end

@testset "JopTranslation, Correctness half width in x" begin
    n1,n2 = 100,120

    m = ones(n1)*[1:n2;]'
    A = JopTranslation(0.0*ones(n1,n2), 0.5*ones(n1,n2))
    d = A*m

    @test m[2:end-1,2:end-1] .+ 0.5 ≈ d[2:end-1,2:end-1]
end
