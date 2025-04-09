using LinearAlgebra, Jets, JetPack, Printf, Test

n1,n2 = 20,5

@testset "JopNormalize, correctness T=$(T), n3=$n3, mode=$mode" for T in (Float64,Float32), n3 in (1,7), mode in (:trace,:shot)
    ϵ = T(1.e-8)
    F = JopNormalize(JetSpace(T,n1,n2,n3); ϵ=ϵ, mode=mode)
    x1 = rand(domain(F))
    x2 = zeros(domain(F))
    Fx1 = F * x1

    for k3 ∈ 1:n3
        (mode === :shot) && (_norm = norm(x1[:,:,k3],2))
        for k2 ∈ 1:n2
            (mode === :trace) && (_norm = norm(x1[:,k2,k3],2))
            x2[:,k2,k3] = x1[:,k2,k3] ./ (_norm + ϵ)
        end
    end
    @test F*x1 ≈ x2
end

@testset "JopNormalize, linearity test, T=$(T), n3=$n3, mode=$mode" for T in (Float32, Float64), n3 in (1,7), mode in (:trace,:shot)
    ϵ = T(1.e-8)
    F = JopNormalize(JetSpace(T,n1,n2,n3);ϵ,mode)
    J = jacobian(F, rand(domain(F)))
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
    lhs,rhs = linearity_test(J')
    @test lhs ≈ rhs
end

@testset "JopNormalize, dot product test, T=$(T), n3=$n3, mode=$mode" for T in (Float32, Float64), n3 in (1,7), mode in (:trace,:shot)
    F = JopNormalize(JetSpace(T,n1,n2,n3))
    J  = jacobian!(F, rand(domain(F)))
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

@testset "JopNormalize, linearization test, T=$(T), n3=$n3, mode=$mode" for T in (Float32, Float64), n3 in (1,7), mode in (:trace,:shot)
    F = JopNormalize(JetSpace(T,n1,n2,n3))
    m0 = rand(domain(F))
    δm = T(0.1) .* rand(domain(F))
    μ = ones(T, 10) 
    μ[1] = 1
    for k = 2:length(μ)
        μ[k] = μ[k-1] / 1.5
    end
    observed, expected = linearization_test(F, m0, μ = μ, δm = δm)
    δ = minimum(abs, observed - expected)
    @test δ < 0.05
end
