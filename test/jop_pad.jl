using JetPack, Jets, Test

@testset "pad, 1D, end" begin
    A = JopPad(JetSpace(Float64, 5), 1:10)
    @test A*[1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test A' ∘ A * [1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 5.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "pad, 1D, start" begin
    A = JopPad(JetSpace(Float64, 5), -4:5)
    @test A*[1.0:5.0;] ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    @test A'∘A*[1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 5.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "pad, 1D, start+end" begin
    A = JopPad(JetSpace(Float64, 5), -4:11)
    @test A*[1.0:5.0;] ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test A'∘A*[1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 5.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "pad, 2D, all sides" begin
    A = JopPad(JetSpace(Float64, 3, 2), 0:5, -1:3)
    @test A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈
    [0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 2.0 0.0;
     0.0 0.0 3.0 4.0 0.0;
     0.0 0.0 5.0 6.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0]
    @test A'∘A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈ [1.0 2.0;3.0 4.0;5.0 6.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "truncate, 2D, all sides" begin
    A = JopPad(JetSpace(Float64, 6, 5), 2:4, 3:4)
    @test A*
    [0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 2.0 0.0;
     0.0 0.0 3.0 4.0 0.0;
     0.0 0.0 5.0 6.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0] ≈
    [1.0 2.0;
     3.0 4.0;
     5.0 6.0]
    @test A'∘A*
    [0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 2.0 0.0;
     0.0 0.0 3.0 4.0 0.0;
     0.0 0.0 5.0 6.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0] ≈
    [0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 2.0 0.0;
     0.0 0.0 3.0 4.0 0.0;
     0.0 0.0 5.0 6.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "pad, 2D, top and sides" begin
    A = JopPad(JetSpace(Float64, 3, 2), -2:3, -1:3)
    @test A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈
    [0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 2.0 0.0;
     0.0 0.0 3.0 4.0 0.0;
     0.0 0.0 5.0 6.0 0.0]
    @test A'∘A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈ [1.0 2.0;3.0 4.0;5.0 6.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "extend pad, 1D, end" begin
    A = JopPad(JetSpace(Float64, 5), 1:10, extend=true)
    @test A*[1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    @test A'∘A*[1.0:5.0;] ≈ [1.0, 2.0, 3.0, 4.0, 30.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "extend pad, 1D, start" begin
    A = JopPad(JetSpace(Float64, 5), -4:5, extend=true)
    @test A*[1.0:5.0;] ≈ [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    @test A'∘A*[1.0:5.0;] ≈ [6.0, 2.0, 3.0, 4.0, 5.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "extend pad, 1D, end" begin
    A = JopPad(JetSpace(Float64, 5), -4:11, extend=true)
    @test A*collect(1.0:5.0) ≈ [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    @test A'∘A*[1.0:5.0;] ≈ [6.0, 2.0, 3.0, 4.0, 35.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "extend pad, 2D" begin
    A = JopPad(JetSpace(Float64, 3, 2), 0:5, -1:3, extend=true)
    @test A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈
    [1.0 1.0 1.0 2.0 2.0;
     1.0 1.0 1.0 2.0 2.0;
     3.0 3.0 3.0 4.0 4.0;
     5.0 5.0 5.0 6.0 6.0;
     5.0 5.0 5.0 6.0 6.0;
     5.0 5.0 5.0 6.0 6.0]
    @test A'∘A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈ [6.0 8.0;9.0 8.0;45.0 36.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end

@testset "extend pad and truc, 2D" begin
    A = JopPad(JetSpace(Float64, 3, 2), -2:3, -1:3, extend=true)
    @test A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈
    [1.0 1.0 1.0 2.0 2.0;
     1.0 1.0 1.0 2.0 2.0;
     1.0 1.0 1.0 2.0 2.0;
     1.0 1.0 1.0 2.0 2.0;
     3.0 3.0 3.0 4.0 4.0;
     5.0 5.0 5.0 6.0 6.0]
    @test A'∘A*[1.0 2.0;3.0 4.0;5.0 6.0] ≈ [12.0 16.0;9.0 8.0;15.0 12.0]
    m = rand(domain(A))
    d = rand(range(A))
    lhs, rhs = dot_product_test(A, m, d)
    @test lhs ≈ rhs
end
