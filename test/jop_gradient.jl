using LinearAlgebra, Jot, Printf, Test

n1,n2,n3 = 11,22,33

@testset "JopGradient2D, linearity test, T=$T" for T in (Float32,Float64)
    A = JopGradient(JetSpace(T,n1,n2), (rand(T),rand(T)))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "Gradient3D, linearity test, T=$(T)" for T in (Float32,Float64)
    A = JopGradient(JetSpace(T,n1,n2,n3), (rand(T),rand(T),rand(T)))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "Gradient2D, dot product test, T=$(T)" for T in (Float32,Float64)
    A = JopGradient(JetSpace(T,n1,n2), (rand(T),rand(T)))
    lhs,rhs = dot_product_test(op, rand(domain(op)), rand(range(op)))
    @test lhs ≈ rhs rtol=1e-6
end

@testset "Gradient3D dot product test, T=$(T)" for T in (Float64,Float32)
    A = JopGradient(JetSpace(T,n1,n2,n3), (rand(T),rand(T),rand(T)))
    lhs,rhs = dot_product_test(A, rand(domain(op)), rand(range(op)))
    @test lhs ≈ rhs rtol=1e-6
end

@testset "Gradient2D correctness test, T=$(T)" for T in (Float64,Float32)
    dom = JetSpace(T,n1,n2)
    delta = (rand(T), rand(T))
    op = JopGradient(dom, delta)
    x = zeros(domain(op))
    α = rand(T)
    β = rand(T)
    for k1 = 1:size(x,1)
        for k2 = 1:size(x,2)
            x[k1,k2] = α * k1 + β * k2
        end
    end
    y = op * x
    diff1 = 0.0
    diff2 = 0.0
    for k1 = 1:size(x,1)
        for k2 = 1:size(x,2)
            if k1 < n1
                diff1 += α - y[k1,k2,1]
            end
            if k2 < n2
                diff2 += β - y[k1,k2,2]
            end
        end
    end
    diff1 /= length(x)
    diff2 /= length(x)
    @test diff1 < 1000 * eps(T)
    @test diff2 < 1000 * eps(T)
end

@testset "Gradient3D correctness test, T=$(T)" for T in (Float64,Float32)
    dom = JetSpace(T,n1,n2,n3)
    delta = (rand(T), rand(T), rand(T))
    op = JopGradient(dom, delta)
    x = zeros(domain(op))
    α = rand(T)
    β = rand(T)
    γ = rand(T)
    for k1 = 1:size(x,1)
        for k2 = 1:size(x,2)
            for k3 = 1:size(x,3)
                x[k1,k2,k3] = α * k1 + β * k2 + γ * k3
            end
        end
    end
    y = op * x
    diff1 = 0.0
    diff2 = 0.0
    diff3 = 0.0
    for k1 = 1:size(x,1)
        for k2 = 1:size(x,2)
            for k3 = 1:size(x,3)
                if k1 < n1
                    diff1 += α - y[k1,k2,k3,1]
                end
                if k2 < n2
                    diff2 += β - y[k1,k2,k3,2]
                end
                if k3 < n3
                    diff3 += γ - y[k1,k2,k3,3]
                end
            end
        end
    end
    diff1 /= length(x)
    diff2 /= length(x)
    diff3 /= length(x)
    @test diff1 < 1000 * eps(T)
    @test diff2 < 1000 * eps(T)
end
