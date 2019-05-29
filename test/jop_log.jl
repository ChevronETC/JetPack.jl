using DSP.Util, LinearAlgebra, Jets, JetPack, Printf, Test

n1,n2 = 33,44

@testset "JopLog" begin

    @testset "JopLog correctness T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
        op = JopLog(JetSpace(T,n1,n2))
        x1 = rand(domain(op)) .+ T(0.0001)
        y1 = op * x1
        y2 = log.(x1)
        dy = y1 .- y2
        diff = sqrt(norm(dy)^2 / length(dy))
        write(stdout, @sprintf("Log calculation -- type(%s) -- rms; %+14.8e\n", T, diff))
        @test diff < 1000 * eps(real(T))
    end

    @testset "Log linearop is linear, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
        op = JopLog(JetSpace(T,n1,n2))
        x = rand(domain(op)) .+ T(0.0001)
        J = jacobian(op, x)

        #  check: J m1 + J m2 - J m3 = J (m1 + m2 - m3)
        m1 = -1 .+ 2 .* rand(domain(op))
        m2 = -1 .+ 2 .* rand(domain(op))
        m3 = -1 .+ 2 .* rand(domain(op))
        dm = (J*m1 .+ J*m2 .- J*m3) .- (J*(m1 .+ m2 .- m3))
        rms = sqrt(norm(dm)^2 / length(dm))
        write(stdout, @sprintf("JotOpLnLog forward is linear -- type(%s) -- rms; %+14.8e\n", T, rms))
        @test rms < 1000 * eps(real(T))

        #  check: J' d1 + J' d2 -J' d3 = J' (d1 + d2 - d3)
        d1 = -1 .+ 2 .* rand(range(op))
        d2 = -1 .+ 2 .* rand(range(op))
        d3 = -1 .+ 2 .* rand(range(op))
        dm = (J'*d1 .+ J'*d2 .- J'*d3) .- (J'*(d1 .+ d2 .- d3))
        rms = sqrt(norm(dm)^2 / length(dm))
        write(stdout, @sprintf("JotOpLnLog adjoint is linear -- type(%s) -- rms; %+14.8e\n", T, rms))
        @test rms < 1000 * eps(real(T))
    end

    @testset "Log dot product test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
        op = JopLog(JetSpace(T,n1,n2))
        x = rand(domain(op)) .+ T(0.0001)
        J  = jacobian(op, x)
        lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(op)), -1 .+ 2 .* rand(range(op)))
        diff = abs((lhs - rhs) / (lhs + rhs))
        write(stdout, @sprintf("Jacobian dot product test -- type(%s) -- rms; %+12.6e\n", T, diff))
        @test diff < 1000 * eps(real(T))
    end

    # note the key here is to increase the size of the nonlinear vector
    @testset "Log linearization test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
        op = JopLog(JetSpace(T,n1,n2))
        m0 = 1000 .* rand(domain(op)) .+ T(0.0001)
        μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096,1/8192])
        dm = rand(domain(op)) .+ T(0.0001)
        observed, expected = linearization_test(op, m0, μ = μ, δm = dm)
        δ = minimum(abs, observed - expected)
        write(stdout, @sprintf("\nLinearization test -- type(%s)\n", T))
        for i = 1:length(observed)
            #write(stdout, @sprintf("mu,observed,expected,diff; %12.6f %12.6f %12.6f %12.6f\n", mu[i], observed[i], expected[i], abs(observed[i] - expected[i])))
        end
        write(stdout, @sprintf("minimum difference %12.6e\n", minimum(abs,observed .- expected)))
        @test δ < 0.1
    end
    nothing
end
nothing
