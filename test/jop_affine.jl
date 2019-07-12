using LinearAlgebra, Jets, JetPack, Test

n1,n2 = 33,44

@testset "JopAffine, correctness T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    b = rand(T,1)
    F = JopAffine(JetSpace(T,n1,n2),b)
    x1 = rand(domain(F))
    @test F*x1 == x1 .+ b

    b1 = rand(T,n1,n2)
    F = JopAffine(b1)
    @test size(domain(F)) == (n1,n2)
    @test size(range(F)) == (n1,n2)
 
    @test F*x1 == x1 .+ b1
end

@testset "JopAffine, linearity test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    @info "Testing scalar b"
    b = rand(T,1)
    F = JopAffine(JetSpace(T,n1,n2),b)
    J = jacobian(F, rand(domain(F)))
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs

    @info "Testing vector b"
    b = rand(T,n1,n2)
    F = JopAffine(b)
    J = jacobian(F, rand(domain(F)))
    lhs,rhs = linearity_test(J)
    @test lhs ≈ rhs
end

@testset "JopAffine, dot product test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    @info "Testing scalar b"
    b = rand(T,1)
    F = JopAffine(JetSpace(T,n1,n2),b)
    J  = jacobian!(F, rand(domain(F)))
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs

    @info "Testing vector b"
    b = rand(T,n1,n2)
    F = JopAffine(b)
     J  = jacobian!(F, rand(domain(F)))
    lhs, rhs = dot_product_test(J, -1 .+ 2 .* rand(domain(J)), -1 .+ 2 .* rand(range(J)))
    @test lhs ≈ rhs
end

#Linearization test modified since the linearization of an affine operator is exactly linear. 
#Usual linearization test will fail since ||F(m + dm) - F(m) - Jdm|| will be close to zero and observed quantity will have 0/0
@testset "JopAffine, linearization test, T=$(T)" for T in (Float64,Float32,Complex{Float64},Complex{Float32})
    @info "Testing scalar b"
    b = rand(T,1)
    F = JopAffine(JetSpace(T,n1,n2),b)
    m0 = rand(domain(F))
    μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096,1/8192])
    δm = rand(domain(F))
    J = jacobian!(F,m0)
    phi = ones(size(μ))
    for i=1:length(μ)
        dnon = F * (m0 .+ μ[i].*δm)
        dlin = F*m0 .+ (μ[i] .* (J*δm))
        phi[i] = norm(dnon .- dlin)
    end
    @show phi
    @test minimum(phi) <= 100*eps(eltype(abs.(b[1])))

    @info "Testing vector b"
    b = rand(T,n1,n2)
    F = JopAffine(b)
    m0 = rand(domain(F))
    μ  = sqrt.([1/1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096,1/8192])
    δm = rand(domain(F))
    J = jacobian!(F,m0)
    phi = ones(size(μ))
    for i=1:length(μ)
        dnon = F * (m0 .+ μ[i].*δm)
        dlin = F*m0 .+ (μ[i] .* (J*δm))
        phi[i] = norm(dnon .- dlin)
    end
    @show phi
    @test minimum(phi) <= 100*eps(eltype(abs.(b[1])))
end
