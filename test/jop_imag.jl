using Jets, JetPack, LinearAlgebra, Printf, Test

n1,n2 = 33,44

@testset "JopImag" begin

	# note: cannot use Jot dot_product_test because we test
	#	imag(domain sum) ≈ range sum
    @testset "ImagPart dot product test, T=$(T)" for T in (Float64,Float32)
        dom = JetSpace(Complex{T},n1,n2)
        op = JopImag(dom)
        lhs,rhs = dot_product_test(op, rand(dom), rand(range(op)))
		@test lhs ≈ rhs
    end

    @testset "ImagPart linearop is linear, T=$(T)" for T in (Float64,Float32)
        dom = JetSpace(Complex{T},n1,n2)
        op = JopImag(dom)
        m1 = -1 .+ 2 .* rand(domain(op))
        m2 = -1 .+ 2 .* rand(domain(op))
        a = op*m1 .+ op*m2
        b = op*(m1 .+ m2)
        diff = sqrt(dot(a-b,a-b) / length(a))
        #write(stdout, @sprintf("linearop is linear test -- type(%s) -- forward rms; %+12.6e\n", T, diff))
        @test diff < 1000 * eps(T)

        d1 = -1 .+ 2 .* rand(range(op))
        d2 = -1 .+ 2 .* rand(range(op))
        a = op'*d1 .+ op'*d2
        b = op'*(d1 .+ d2)
        diff = sqrt(abs(dot(a-b,a-b)) / length(a))
        #write(stdout, @sprintf("linearop is linear test -- type(%s) -- adjoint rms; %+12.6e\n", T, diff))
        @test diff < 10 * eps(T)
    end

    @test_skip @testset "ImagPart correctness test, T=$(T)" for T in (Float64,Float32)
        dom = JetSpace(Complex{T},n1,n2)
        rng = JetSpace(T,n1,n2)
        op = JotOpImagPart(dom,rng)
        x = rand(domain(op))
        a = op * x
		b = 1/(2*im) .* (x .- conj.(x))
        diff = sqrt(abs(dot(a-b,a-b)) / length(a))
        #write(stdout, @sprintf("correctness test -- type(%s) -- diff; %+12.6e\n", T, diff))
        @test diff < 10 * eps(T)
    end
    nothing
end
nothing
