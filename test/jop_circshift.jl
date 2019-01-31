using JetPack, Jets, Test

n1,n2,n3 = 10,11,12

@testset "JopCircShift" begin
    @testset "CircShift linearop is linear, size=$(N), T=$(T)" for T in (Float64,Float32), N in ((n1,n2), (n1,n2,n3))
        op = JopCircShift(JetSpace(T,N), (1,2,3)[1:length(N)])
        m1 = -1 .+ 2 .* rand(domain(op))
        m2 = -1 .+ 2 .* rand(domain(op))
        @test (op*m1 .+ op*m2) ≈ op*(m1 .+ m2)
        d1 = -1 .+ 2 .* rand(range(op))
        d2 = -1 .+ 2 .* rand(range(op))
        @test (op'*d1 .+ op'*d2) ≈ op'*(d1 .+ d2)
    end

    @testset "CircShift dot product test, size=$(N), T=$(T)" for T in (Float64,Float32), N in ((n1,n2), (n1,n2,n3))
        op = JopCircShift(JetSpace(T,N), (1,2,3)[1:length(N)])
        lhs,rhs = dot_product_test(op, rand(domain(op)), rand(range(op)))
        @test isapprox(lhs, rhs, rtol=1e-6)
    end
end
