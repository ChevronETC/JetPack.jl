using JetPack, Jets, Test

n1dom2d,n2dom2d = 10,11
n1rng2d,n2rng2d = 29,31

@testset "JopInterp 2D, interpolation test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom2d,n2dom2d), JetSpace(T,n1rng2d,n2rng2d))
    xdom = zeros(domain(A))
    a,b = rand(2)
    for k1 ∈ 1:n1dom2d
        x1 = (k1-1) / (n1dom2d-1)
        for k2 ∈ 1:n2dom2d
            x2 = (k2-1) / (n2dom2d-1)
            xdom[k1,k2] = a*x1 + b*x2
        end
    end
    xrng = A * xdom
    for k1 ∈ 1:n1rng2d
        x1 = (k1-1) / (n1rng2d-1)
        for k2 ∈ 1:n2rng2d
            x2 = (k2-1) / (n2rng2d-1)
            @test xrng[k1,k2] ≈ a*x1 + b*x2
        end
    end
end

@testset "JopInterp 2D, linearity test" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom2d,n2dom2d), JetSpace(T,n1rng2d,n2rng2d))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopInterp 2D, dot product test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom2d,n2dom2d), JetSpace(T,n1rng2d,n2rng2d))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopInterp 2D, test bugfix for specific dom,rng, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,280,222), JetSpace(T,281,2221))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

n1dom3d,n2dom3d,n3dom3d = 5,6,7
n1rng3d,n2rng3d,n3rng3d = 13,15,17

@testset "JopInterp 3D, interpolation test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom3d,n2dom3d,n3dom3d), JetSpace(T,n1rng3d,n2rng3d,n3rng3d))
    xdom = zeros(domain(A))
    a,b,c = rand(3)
    for k1 ∈ 1:n1dom3d
        x1 = (k1-1) / (n1dom3d-1)
        for k2 ∈ 1:n2dom3d
            x2 = (k2-1) / (n2dom3d-1)
            for k3 ∈ 1:n3dom3d
                x3 = (k3-1) / (n3dom3d-1)
                xdom[k1,k2,k3] = a*x1 + b*x2 + c*x3
            end
        end
    end
    xrng = A * xdom
    for k1 ∈ 1:n1rng3d
        x1 = (k1-1) / (n1rng3d-1)
        for k2 ∈ 1:n2rng3d
            x2 = (k2-1) / (n2rng3d-1)
            for k3 ∈ 1:n3rng3d
                x3 = (k3-1) / (n3rng3d-1)
                @test xrng[k1,k2,k3] ≈ a*x1 + b*x2 + c*x3
            end
        end
    end
end

@testset "JopInterp 3D, linearity test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom3d,n2dom3d,n3dom3d), JetSpace(T,n1rng3d,n2rng3d,n3rng3d))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopInterp 3D, dot product test, T=$(T)" for T in (Float32,Float64)
    A = JopInterp(JetSpace(T,n1dom3d,n2dom3d,n3dom3d), JetSpace(T,n1rng3d,n2rng3d,n3rng3d))
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
