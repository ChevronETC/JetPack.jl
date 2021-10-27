using JetPack, Jets, Test

@testset "JopZero, zero test suite" begin
    for T1 in (Float32,Float64,ComplexF32,ComplexF64), T2 in (Float32,Float64,ComplexF32,ComplexF64), n1 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5)), n2 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5))
        A = JopZero(JetSpace(T1,n1...), JetSpace(T2,n2...))
        @test zeros(JetSpace(T2,n2...)) ≈ A*rand(JetSpace(T1,n1...))
        @test zeros(JetSpace(T1,n1...)) ≈ A'*rand(JetSpace(T2,n2...))
        close(A)
    end
end

@testset "JopZero, linearity test suite" begin
    for T1 in (Float32,Float64,ComplexF32,ComplexF64), T2 in (Float32,Float64,ComplexF32,ComplexF64), n1 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5)), n2 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5))
        A = JopZero(JetSpace(T1,n1...), JetSpace(T2,n2...))
        lhs,rhs = linearity_test(A)
        @test lhs ≈ rhs
        close(A)
    end
end

@testset "JopZero, dot product test suite" begin
    for T1 in (Float32,Float64,ComplexF32,ComplexF64), T2 in (Float32,Float64,ComplexF32,ComplexF64), n1 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5)), n2 in ((4,), (1,2), (2,5,1), (23,2,1,4), (6,2,3,1,3,5))
        A = JopZero(JetSpace(T1,n1...), JetSpace(T2,n2...))
        lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
        @test lhs ≈ rhs
        close(A)
    end
end

nothing
