"""
    A = JopRBF(dom, rng, nodes; delta)

Meshless interpolation from a *scattered* set of control nodes to a regular fine
grid, using **normalized compact-support radial basis functions** (Wendland C2).

This is a reduced (coarse) model parameterization: the domain coefficients live on
an arbitrary cloud of points (not a grid), and the operator evaluates a smooth
field from them on the range grid. It is built for very large fine grids (e.g.
`1000^3`): every fine point depends on only the O(1) nodes within a support
radius, so the forward and adjoint are applied matrix-free with a node bucket
index; nothing dense is formed or stored.

The field is the normalized (Shepard) combination
```
         Σ_j φ(r_j(x)) c_j
d(x) = ─────────────────────,   r_j(x) = sqrt( Σ_k ((x_k - ξ_{j,k}) / δ_{k,j})^2 )
          Σ_j φ(r_j(x))
```
with the Wendland C2 kernel `φ(r) = (1-r)^4 (4r+1)` on `0 ≤ r < 1` (zero beyond).
Properties, all exercised by the unit tests:

* **Compact support / scalable** -- `φ` vanishes for `r ≥ 1`, so node `j` only
  influences fine points inside the ellipsoid of radii `δ_{:,j}` about `ξ_j`.
* **C2 smooth** -- `φ`, `φ'` and `φ''` vanish at `r = 1`, so the field is twice
  continuously differentiable (no bilinear-style creases / hot zones).
* **Reproduces constants** -- the normalization is a partition of unity, so a
  constant coefficient vector maps to that constant wherever there is coverage.
* **No overshoot** -- `φ ≥ 0`, so the field is a convex combination of nearby
  coefficients and stays within their min/max (no ringing / halos).

# Arguments

`dom`: operator domain, `JetSpace(T, M)`, the `M` scattered node coefficients.\\
`rng`: operator range, `JetSpace(T, n1[, n2[, n3]])`, the fine grid. The number of
active spatial dimensions is `ndims(rng)` (1, 2 or 3).\\
`nodes`: a `ndims(rng) x M` matrix of node coordinates in *fine-grid index
coordinates* (fine sample `i` along an axis sits at coordinate `i`, so
coordinates lie in `[1, n_k]`). Placing some nodes on the boundary of the model
guarantees every fine point is covered (no extrapolation holes).\\
`delta`: the compact-support radius, in fine-grid index units. It may be

* a **scalar** (isotropic, the same radius for every node and axis);
* a length-`ndims(rng)` **vector** (per axis, e.g. for anisotropic grid spacing);
  or
* a `ndims(rng) x M` **matrix** (per node and axis). This is the multiresolution
  form: each node `j` carries its own support radii `δ_{:,j}`, so a
  points-per-wavelength node cloud (dense with small support near the surface,
  sparse with large support at depth) resolves fine detail where nodes are dense
  and stays smooth where they are sparse.

A fine point is influenced by node `j` when it is within `δ_{:,j}` of `ξ_j` along
every axis. The support must be large enough that every fine point of interest is
within some node's support; fine points not covered by any node (`Σ_j φ = 0`) map
to zero.

# Examples

```julia
# 2D: 300 scattered nodes onto a 400 x 500 grid
nz, nx, M = 400, 500, 300
nodes = vcat((1 .+ (nz-1).*rand(1,M)), (1 .+ (nx-1).*rand(1,M)))   # 2 x M, index coords
A = JopRBF(JetSpace(Float64,M), JetSpace(Float64,nz,nx), nodes; delta = 60.0)
c = rand(domain(A))
d = A*c
```

```julia
# per-node (multiresolution) support: each node gets its own radii (2 x M)
deltas = vcat(fill(15.0, 1, M), fill(15.0, 1, M))   # here isotropic, but may vary per node
A = JopRBF(JetSpace(Float64,M), JetSpace(Float64,nz,nx), nodes; delta = deltas)
```
"""
function JopRBF(dom::JetSpace{T,1}, rng::JetSpace{T,D}, nodes::AbstractMatrix; delta) where {T,D}
    (D >= 1 && D <= 3) || error("JopRBF supports 1, 2 or 3 spatial (range) dimensions, got $D")
    M = size(dom)[1]
    n = size(rng)
    size(nodes) == (D, M) ||
        error("nodes must be (ndims(range)=$D) x (length(domain)=$M), got $(size(nodes))")

    δ = _delta_matrix(delta, D, M)                     # (D, M) per-node per-axis radii
    all(x -> x > 0, δ) || error("delta must be positive")

    # bucket size = the largest support radius per axis, so that any node covering
    # a fine point sits within one bucket of it (exact neighbor search).
    bucketsize = ntuple(k -> maximum(@view δ[k, :]), D)

    nodesf = Matrix{Float64}(nodes)
    buckets = _build_buckets(nodesf, bucketsize)
    invw = _build_invw(n, nodesf, δ, bucketsize, buckets)

    JopLn(; dom, rng, df! = JopRBF_df!, df′! = JopRBF_df′!,
        s = (; nodes = nodesf, delta = δ, bucketsize, buckets, invw))
end

export JopRBF

_delta_matrix(delta::Number, D, M) = fill(Float64(delta), D, M)
function _delta_matrix(delta::AbstractVector, D, M)
    length(delta) == D || error("delta vector length ($(length(delta))) must equal number of range dimensions ($D)")
    Float64[delta[k] for k = 1:D, j = 1:M]
end
function _delta_matrix(delta::AbstractMatrix, D, M)
    size(delta) == (D, M) || error("delta matrix must be (ndims(range)=$D) x (length(domain)=$M), got $(size(delta))")
    Matrix{Float64}(delta)
end

"""
Wendland C2 radial basis function: `φ(r) = (1-r)^4 (4r+1)` for `0 ≤ r < 1`, else
`0`. Positive definite in up to three dimensions, non-negative, and C2 (value,
first and second derivatives vanish at `r = 1`).
"""
@inline function wendland_c2(r::Float64)
    if r < 1.0
        om = 1.0 - r
        return om * om * om * om * (4.0 * r + 1.0)
    else
        return 0.0
    end
end

# assign each node to a bucket of size `bucketsize` (per axis) for O(1) neighbor search
function _build_buckets(nodes::Matrix{Float64}, bucketsize::NTuple{D,Float64}) where {D}
    M = size(nodes, 2)
    buckets = Dict{NTuple{D,Int}, Vector{Int}}()
    for j = 1:M
        b = ntuple(k -> floor(Int, nodes[k, j] / bucketsize[k]), D)
        push!(get!(() -> Int[], buckets, b), j)
    end
    buckets
end

# precompute the normalization: invw[x] = 1 / Σ_j φ(r_j(x)) (0 where uncovered).
# Built once from node positions and reused for every apply.
function _build_invw(n::NTuple{D,Int}, nodes::Matrix{Float64}, delta::Matrix{Float64}, bucketsize::NTuple{D,Float64}, buckets) where {D}
    invw = zeros(Float64, n)
    Threads.@threads :static for I in CartesianIndices(invw)
        x = Tuple(I)
        w = _accumulate_weight(x, nodes, delta, bucketsize, buckets)
        @inbounds invw[I] = w > 0.0 ? 1.0 / w : 0.0
    end
    invw
end

# sum of kernel weights at fine coordinate `x` over nodes in the 3^D neighbor buckets
@inline function _accumulate_weight(x::NTuple{D,Int}, nodes::Matrix{Float64}, delta::Matrix{Float64}, bucketsize::NTuple{D,Float64}, buckets) where {D}
    bx = ntuple(k -> floor(Int, x[k] / bucketsize[k]), D)
    w = 0.0
    for off in CartesianIndices(ntuple(_ -> -1:1, D))
        ob = Tuple(off)
        b = ntuple(k -> bx[k] + ob[k], D)
        nb = get(buckets, b, nothing)
        nb === nothing && continue
        for j in nb
            r2 = 0.0
            @inbounds for k = 1:D
                dk = (x[k] - nodes[k, j]) / delta[k, j]
                r2 += dk * dk
            end
            if r2 < 1.0
                w += wendland_c2(sqrt(r2))
            end
        end
    end
    w
end

#
# forward: thread over fine points, gather nearby nodes via the bucket index
# (each fine point writes only its own output, so no data races).
#
function JopRBF_df!(d::AbstractArray{T,D}, c::AbstractVector{T}; nodes, delta, bucketsize, buckets, invw, kwargs...) where {T,D}
    Threads.@threads :static for I in CartesianIndices(d)
        x = Tuple(I)
        bx = ntuple(k -> floor(Int, x[k] / bucketsize[k]), D)
        num = 0.0
        for off in CartesianIndices(ntuple(_ -> -1:1, D))
            ob = Tuple(off)
            b = ntuple(k -> bx[k] + ob[k], D)
            nb = get(buckets, b, nothing)
            nb === nothing && continue
            for j in nb
                r2 = 0.0
                @inbounds for k = 1:D
                    dk = (x[k] - nodes[k, j]) / delta[k, j]
                    r2 += dk * dk
                end
                if r2 < 1.0
                    @inbounds num += wendland_c2(sqrt(r2)) * c[j]
                end
            end
        end
        @inbounds d[I] = convert(T, num * invw[I])
    end
    d
end

#
# adjoint: thread over nodes, iterate each node's fine-grid bounding box (each
# node writes only its own coefficient, so no data races). This enumerates the
# same (fine point, node) pairs as the forward, so it is the exact transpose.
#
function JopRBF_df′!(c::AbstractVector{T}, d::AbstractArray{T,D}; nodes, delta, bucketsize, buckets, invw, kwargs...) where {T,D}
    n = size(d)
    Threads.@threads :static for j = 1:length(c)
        xlo = ntuple(k -> max(1, ceil(Int, nodes[k, j] - delta[k, j])), D)
        xhi = ntuple(k -> min(n[k], floor(Int, nodes[k, j] + delta[k, j])), D)
        acc = 0.0
        for I in CartesianIndices(ntuple(k -> xlo[k]:xhi[k], D))
            x = Tuple(I)
            r2 = 0.0
            @inbounds for k = 1:D
                dk = (x[k] - nodes[k, j]) / delta[k, j]
                r2 += dk * dk
            end
            if r2 < 1.0
                @inbounds acc += wendland_c2(sqrt(r2)) * (convert(Float64, d[I]) * invw[I])
            end
        end
        @inbounds c[j] = convert(T, acc)
    end
    c
end
