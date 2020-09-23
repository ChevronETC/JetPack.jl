using JetPack, Jets, Test, FFTW

n1,n2 = 101,201

@testset "JopRemoveDC, correctness test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64}) 
    A = JopRemoveDC(JetSpace(T,n1,n2))
    x = -1 .+ 2 .* rand(domain(A)) .+ rand(T)
    y2 = A * x
    X = fft(x,(1,))
    X[1,:] .= 0
    y1 = ifft(X,(1,))
    @test y1 ≈ y2
end

@testset "JopRemoveDC, linearity test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64}) 
    A = JopRemoveDC(JetSpace(T,n1,n2))
    lhs,rhs = linearity_test(A)
    @test lhs ≈ rhs
end

@testset "JopRemoveDC, dot product test, T=$T" for T in (Float32,Float64,Complex{Float32},Complex{Float64})
    #=
    I observe the following bug on windows:
    x = rand(Complex{Float64},101,201)
    using LinearAlgebra
    dot(x,x) # crashes Julia, presumably due to an issue with BLAS on windows.
    =#
    if Sys.iswindows() && T == Complex{Float64}
        @warn "dot product test fails on Windows due to the LinearAlgebra.dot method not working for Complex{Float64}"
    else
        A = JopRemoveDC(JetSpace(T,n1,n2))
        lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
        @test lhs ≈ rhs
    end
end
