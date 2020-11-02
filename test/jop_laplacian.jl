using JetPack, Jets, Test

@testset "JopLaplacian, dot product T=$(T)" for T in (Float32, Float64)
    sp = JetSpace(T,64,128)
    A = JopLaplacian(sp)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs

    sp = JetSpace(T,16,32,64)
    A = JopLaplacian(sp)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

# note this is a second order derivative, so we loosen up tolerance from 
# 2nd to 4th root epsilon(type)
@testset "JopLaplacian, action T=$(T)" for T in (Float32, Float64)
    tol = eps(T)^0.25

    sp = JetSpace(T,64,128)
    A = JopLaplacian(sp)
    a,b = rand(T),rand(T)
    x = zeros(sp)
    nz,nx=size(x)
    for ix=1:nx, iz=1:nz
        x[iz,ix] = a * iz^2 + b * ix^2
    end
    y = A * x
    expected = 2 * (a + b)
    z = y[2:nz-1,2:nx-1] .- expected .* ones(T,nz-2,nx-2)
    @test maximum(abs,z) < tol

    sp = JetSpace(T,16,32,64)
    A = JopLaplacian(sp)
    a,b,c = rand(T),rand(T),rand(T)
    x = zeros(sp)
    nz,ny,nx=size(x)
    for ix=1:nx, iy=1:ny, iz=1:nz
        x[iz,iy,ix] = a * iz^2 + b * iy^2 + c * iz^2
    end
    y = A * x
    expected = 2 * (a + b + c)
    z = y[2:nz-1,2:ny-1,2:nx-1] .- expected .* ones(T,nz-2,ny-2,nx-2)
    @test maximum(abs,z) < tol
end