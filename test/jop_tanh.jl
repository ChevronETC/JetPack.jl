using LinearAlgebra, Jets, JetPack, Test

n1,n2 = 33,44

@testset "JopTanh, correctness T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopTanh(JetSpace(T,n1,n2))
    x1 = rand(domain(F)) 
    @test F*x1 ≈ tanh.(x1)
end

@testset "JopTanh, linearity test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopTanh(JetSpace(T,n1,n2))
    J = jacobian(F, rand(domain(F)) )
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
end

@testset "JopTanh, dot product test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopTanh(JetSpace(T,n1,n2))
    J  = jacobian!(F, rand(domain(F)) )
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

# note the key here is to increase the size of the nonlinear vector
@testset "JopTanh, linearization test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopTanh(JetSpace(T,n1,n2))
    m0 = -0.5 .* rand(domain(F))
    μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048])
    δm = rand(domain(F))
    observed, expected = linearization_test(F, m0, μ = μ, δm = δm)
    δ = minimum(abs, observed - expected)
    @test δ < 0.1
end

@testset "JopTanh, odd funcion T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    F = JopTanh(JetSpace(T,n1,n2))
    x1 = 2 .* (-1 .+ 2 .* rand(domain(F)))
    @test F * (-1 .* x1) ≈ -(F * x1)
end
