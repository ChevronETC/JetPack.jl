using JetPack, Jets, Test

@testset "JopDerivative, 2D, along first dimension" begin
    x = range(0, stop=4*pi, length=256)
    m = sin.(x) * ones(1,10)
    A = JopDerivative(JetSpace(Float64,256,10), delta=x[2]-x[1])
    d = A*m
    a = A'*d

    dim=1
    Rpre = CartesianIndices(size(m)[1:dim-1])
    Rpost = CartesianIndices(size(m)[dim+1:end])
    n = size(m, dim)

    d_expected = cos.(x) * ones(1,10)
    for j = 1:10, i = 5:length(x)-5
        @test isapprox(d[i,j], d_expected[i,j], rtol=1e-5)
    end

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopDerivative, 3D, along first dimension" begin
    x = range(0, stop=4*pi, length=256)
    m = zeros(256,10,1)
    m[:,:,:] = sin.(x) * ones(1,10)
    A = JopDerivative(JetSpace(Float64,256,10,1), delta=x[2]-x[1])
    d = A*m
    a = A'*d

    d_expected = cos.(x) * ones(1,10)
    for j = 1:10, i = 5:length(x)-5
        @test isapprox(d[i,j,1], d_expected[i,j,1], rtol=1e-5)
    end

    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end
