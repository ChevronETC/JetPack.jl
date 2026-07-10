using Jets, LinearAlgebra, Test, JetPack

# a regular lattice of nodes (in fine-index coordinates) that includes the domain
# boundary, optionally jittered, guaranteeing full coverage for delta >= 1.5*h.
function lattice_nodes(n::NTuple{D,Int}, h::Real; jitter = 0.0, rng = nothing) where {D}
    axes1d = ntuple(k -> unique(vcat(collect(1.0:h:n[k]), Float64(n[k]))), D)
    grid = Iterators.product(axes1d...)
    pts = [collect(Float64, p) for p in grid]
    nodes = reduce(hcat, pts)
    if jitter > 0
        r = rng === nothing ? () : (rng,)
        for j = 1:size(nodes, 2), k = 1:D
            interior = nodes[k, j] > 1.5 && nodes[k, j] < n[k] - 0.5
            interior && (nodes[k, j] += jitter * (rand(r...) - 0.5))
        end
    end
    nodes
end

@testset "JopRBF, wendland_c2 kernel is C2 and non-negative" begin
    φ = JetPack.wendland_c2
    @test φ(0.0) == 1.0
    @test φ(1.0) == 0.0
    @test φ(1.5) == 0.0
    @test all(φ(r) >= 0 for r in range(0, 1, length = 101))
    # value, first and second derivative all vanish at r -> 1^- (=> C2 at support edge)
    h = 1e-4
    r0 = 1.0 - h
    d1 = (φ(r0 + h/2) - φ(r0 - h/2)) / h
    d2 = (φ(r0 + h) - 2φ(r0) + φ(r0 - h)) / h^2
    @test isapprox(φ(r0), 0.0; atol = 1e-8)
    @test isapprox(d1, 0.0; atol = 1e-3)
    @test isapprox(d2, 0.0; atol = 1e-1)
end

@testset "JopRBF, dot product test - $(D)D" for D in (1, 2, 3)
    n = D == 1 ? (151,) : D == 2 ? (61, 47) : (23, 19, 17)
    M = 40
    nodes = zeros(D, M)
    for k = 1:D, j = 1:M
        nodes[k, j] = 1 + (n[k] - 1) * rand()
    end
    delta = 0.5 * minimum(n)
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = delta)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopRBF, reproduces constants (partition of unity) - $(D)D" for D in (1, 2, 3)
    n = D == 1 ? (120,) : D == 2 ? (60, 50) : (24, 20, 18)
    h = D == 1 ? 12 : D == 2 ? 12 : 6
    nodes = lattice_nodes(n, h)
    M = size(nodes, 2)
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = 1.5h, precondition = false)
    d = A * ones(domain(A))
    @test all(x -> isapprox(x, 1.0; atol = 1e-12), d)   # full coverage => constant reproduced everywhere
end

@testset "JopRBF, no overshoot (convex combination) - $(D)D" for D in (2, 3)
    n = D == 2 ? (60, 50) : (24, 20, 18)
    h = D == 2 ? 12 : 6
    nodes = lattice_nodes(n, h)
    M = size(nodes, 2)
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = 1.5h, precondition = false)
    c = rand(domain(A))                     # coefficients in [0,1]
    d = A * c
    @test all(x -> -1e-12 <= x <= 1 + 1e-12, d)
    @test minimum(d) >= minimum(c) - 1e-12
    @test maximum(d) <= maximum(c) + 1e-12
end

@testset "JopRBF, compact support - influence vanishes beyond delta" begin
    n = (80, 90)
    nodes = [20.0 60.0;    # node 1 at (20,30), node 2 at (60,70)
             30.0 70.0]
    delta = 18.0
    A = JopRBF(JetSpace(Float64, 2), JetSpace(Float64, n...), nodes; delta = delta)
    e1 = zeros(domain(A)); e1[1] = 1.0
    d = A * e1
    for I in CartesianIndices(d)
        x = Tuple(I)
        r = sqrt((x[1] - nodes[1,1])^2 + (x[2] - nodes[2,1])^2) / delta
        if r >= 1.0
            @test d[I] == 0.0                # zero outside node 1's support
        end
    end
    @test any(d .> 0)                        # nonzero somewhere inside the support
end

@testset "JopRBF, anisotropic (per-axis) delta - dot product test" begin
    n = (50, 40)
    M = 30
    nodes = zeros(2, M)
    for j = 1:M
        nodes[1, j] = 1 + (n[1] - 1) * rand()
        nodes[2, j] = 1 + (n[2] - 1) * rand()
    end
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = [10.0, 20.0])
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
end

@testset "JopRBF, per-node (multiresolution) delta" begin
    # each node carries its own support radii (2 x M). Small support near the
    # surface, large support at depth (a points-per-wavelength style cloud).
    n = (120, 90)
    nodes = lattice_nodes(n, 12)
    M = size(nodes, 2)
    # per-node radius grows with depth (row-1 coordinate)
    deltas = zeros(2, M)
    for j = 1:M
        z = nodes[1, j]
        rad = 14.0 + 40.0 * (z - 1) / (n[1] - 1)   # 14 shallow -> 54 deep
        deltas[1, j] = rad
        deltas[2, j] = rad
    end
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = deltas, precondition = false)
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
    # constants still reproduced where covered
    d = A * ones(domain(A))
    covered = d .> 0
    @test all(x -> isapprox(x, 1.0; atol = 1e-12), d[covered])
    @test count(covered) == length(d)            # full coverage for these radii
    # a shallow node influences a smaller region than a deep node
    jshallow = argmin(nodes[1, :]); jdeep = argmax(nodes[1, :])
    supp(j) = (e = zeros(domain(A)); e[j] = 1.0; count(!iszero, A * e))
    @test supp(jdeep) > supp(jshallow)
end

@testset "JopRBF, single precision support" begin
    n = (60, 50)
    nodes = lattice_nodes(n, 12)
    M = size(nodes, 2)
    A = JopRBF(JetSpace(Float32, M), JetSpace(Float32, n...), Float32.(nodes); delta = 18.0f0, precondition = false)
    @test eltype(domain(A)) == Float32
    lhs, rhs = dot_product_test(A, rand(domain(A)), rand(range(A)))
    @test lhs ≈ rhs
    d = A * ones(domain(A))
    @test all(x -> isapprox(x, 1f0; atol = 1f-5), d)
end

@testset "JopRBF, precondition=true equalizes column norms (multiresolution)" begin
    # a node cloud with support that varies ~4x with depth => the raw kernel has a big column-norm (support-
    # size) imbalance; the default precondition=true rescales each coefficient by 1/‖A e_j‖ so the columns are
    # ~equal norm (order unity), while keeping an exact adjoint. This is the Jacobi preconditioner for a
    # reduced-parameterization inverse problem, moved into the operator at construction.
    n = (120, 90)
    nodes = lattice_nodes(n, 12)
    M = size(nodes, 2)
    deltas = zeros(2, M)
    for j = 1:M
        rad = 14.0 + 40.0 * (nodes[1, j] - 1) / (n[1] - 1)   # 14 shallow -> 54 deep
        deltas[1, j] = rad; deltas[2, j] = rad
    end
    colnorms(A) = [norm(A * setindex!(zeros(domain(A)), 1.0, j)) for j in 1:M]
    spread(v) = maximum(v) / minimum(v)
    A0 = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = deltas, precondition = false)
    A1 = JopRBF(JetSpace(Float64, M), JetSpace(Float64, n...), nodes; delta = deltas, precondition = true)
    @test spread(colnorms(A0)) > 3.0              # raw kernel: large support-size imbalance
    @test spread(colnorms(A1)) < 1.5              # preconditioned: ~order unity
    lhs, rhs = dot_product_test(A1, rand(domain(A1)), rand(range(A1)))
    @test lhs ≈ rhs                               # exact adjoint preserved under the preconditioner
end

nothing
