using Jets, LinearAlgebra, IterativeSolvers, Test, JetPack

PLOTS = parse(Int32, get(ENV, "JP_PLOTS", "0"))

n1 = 101
n2 = 60
n3 = 31
n4 = 5

nc1 = 51
nc2 = 30
nc3 = 31

maxiter = 5

@testset "JopCubicSpline, dot product test - N = $(N)" for N = 1:4

    dom, rng = nothing, nothing
    if N == 1
        dom = JetSpace(Float32, nc1)
        rng = JetSpace(Float32, n1)
    elseif N == 2
        dom = JetSpace(Float32, nc1, nc2)
        rng = JetSpace(Float32, n1, n2)
    elseif N == 3
        dom = JetSpace(Float32, nc1, nc2, nc3)
        rng = JetSpace(Float32, n1, n2, n3)
    elseif N == 4
        dom = JetSpace(Float32, nc1, nc2, nc3, n4)
        rng = JetSpace(Float32, n1, n2, n3, n4)
    end
    
    A = JopCubicSpline(dom, rng)
    lhs,rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @show lhs,rhs,(lhs-rhs)/(lhs+rhs)
    @test lhs ≈ rhs
end


@testset "JopCubicSpline, smoothing test - 1D" begin

    dom = JetSpace(Float32, nc1)
    rng = JetSpace(Float32, n1)
    factor = (n1 - 1) / (nc1 - 1)
    
    A = JopCubicSpline(dom, rng)

    x = LinRange(0,1,n1)
    xdec = LinRange(0,1,nc1)
    d = zeros(Float32, n1)
    d .= cos.(8 * π .* x)
    m = convert(Matrix, A) \ d
    dsmth = A * m

    @test norm(d - dsmth) / norm(d) < 0.01

    if PLOTS == 1
        using PyPlot
        close("all")

        figure(figsize=(8,4))
        plot(x, d, label="Original",linewidth=3,color="black")
        plot(xdec, m, label="Spline",linewidth=2,color="blue")
        plot(x, dsmth, label="Smoothed",linewidth=1,color="red")
        xlim([0,1])
        xlabel("x")
        ylabel("y")
        title("JopCubicSpline - 1D")
        legend()
        tight_layout()
        savefig("smoothing1D.JopCubicSpline.png", dpi=100) 
    end
end

@testset "JopCubicSpline, smoothing test - 2D" begin

    dom = JetSpace(Float32, nc1, nc2)
    rng = JetSpace(Float32, n1, n2)
    factor = (n1 - 1) * (n2 - 1) / ((nc1 - 1) * (nc2 - 1))
    
    A = JopCubicSpline(dom, rng)

    x = LinRange(0,1,n1)
    y = LinRange(0,1,n2)
    xdec = LinRange(0,1,nc1)
    ydec = LinRange(0,1,nc2)
    d = zeros(Float32, n1, n2)
    d .= (cos.(8 * π .* x)) * (cos.(8 * π .* y))'
    m = reshape(convert(Matrix, A) \ vec(d), domain(A))
    dsmth = A * m

    @test norm(d - dsmth) / norm(d) < 0.01

    if PLOTS == 1
        using PyPlot
        close("all")

        figure(figsize=(8,8))
        subplot(2,2,1)
        imshow(d, aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Original")
        subplot(2,2,2)
        imshow(m, aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Spline")
        subplot(2,2,3)
        imshow(dsmth, aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Smoothed")
        tight_layout()
        savefig("smoothing2D.JopCubicSpline.png", dpi=100) 
    end
end

@testset "JopCubicSpline, smoothing test - 3D" begin

    dom = JetSpace(Float32, nc1, nc2, nc3)
    rng = JetSpace(Float32, n1, n2, n3)
    factor = (n1 - 1) * (n2 - 1) * (n3 - 1) / ((nc1 - 1) * (nc2 - 1) * (nc3 - 1))
    
    A = JopCubicSpline(dom, rng)

    x = LinRange(0,1,n1)
    y = LinRange(0,1,n2)
    z = LinRange(0,1,n3)
    xdec = LinRange(0,1,nc1)
    ydec = LinRange(0,1,nc2)
    zdec = LinRange(0,1,nc3)
    d = zeros(Float32, n1, n2, n3)
    for iz = 1:n3
        d[:,:,iz] .= (cos.(8 * π .* x)) * (cos.(8 * π .* y))' .* cos.(8 * π .* z[iz])
    end
    m = A' * d
    m ./= factor
    lsqr!(vec(m), vec(A), vec(d); maxiter = maxiter, verbose = true)
    dsmth = A * m
    @test norm(d - dsmth) / norm(d) < 0.01

    if PLOTS == 1
        using PyPlot
        close("all")

        figure(figsize=(8,8))
        subplot(2,2,1)
        imshow(d[:,:,div(n3,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Original")
        subplot(2,2,2)
        imshow(m[:,:,div(nc3,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Spline")
        subplot(2,2,3)
        imshow(dsmth[:,:,div(n3,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Smoothed")
        tight_layout()
        savefig("smoothing3D.JopCubicSpline.png", dpi=100)
    end
end

@testset "JopCubicSpline, smoothing test - 3D with multicomponents" begin

    dom = JetSpace(Float32, nc1, nc2, nc3, n4)
    rng = JetSpace(Float32, n1, n2, n3, n4)
    factor = (n1 - 1) * (n2 - 1) * (n3 - 1) / ((nc1 - 1) * (nc2 - 1) * (nc3 - 1))
    
    A = JopCubicSpline(dom, rng)

    x = LinRange(0,1,n1)
    y = LinRange(0,1,n2)
    z = LinRange(0,1,n3)
    xdec = LinRange(0,1,nc1)
    ydec = LinRange(0,1,nc2)
    zdec = LinRange(0,1,nc3)
    d = zeros(Float32, n1, n2, n3, n4)
    for iz = 1:n3
        d[:,:,iz,1] .= (cos.(8 * π .* x)) * (cos.(8 * π .* y))' .* cos.(8 * π .* z[iz])
    end
    for c = 2:n4
        d[:,:,:,c] .= d[:,:,:,1]
    end
    m = A' * d
    m ./= factor
    lsqr!(vec(m), vec(A), vec(d); maxiter = maxiter, verbose = true)
    dsmth = A * m

    @test norm(d - dsmth) / norm(d) < 0.01

    if PLOTS == 1
        using PyPlot
        close("all")

        figure(figsize=(8,8))
        subplot(2,2,1)
        imshow(d[:,:,div(n3,2),div(n4,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Original")
        subplot(2,2,2)
        imshow(m[:,:,div(nc3,2),div(n4,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Spline")
        subplot(2,2,3)
        imshow(dsmth[:,:,div(n3,2),div(n4,2)], aspect="auto", extent=[0,1,1,0], clim=[-1, 1], cmap="seismic")
        title("Smoothed")
        tight_layout()
        savefig("smoothing3D-MC.JopCubicSpline.png", dpi=100) 
    end
end

nothing