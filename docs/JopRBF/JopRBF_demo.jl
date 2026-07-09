# Kernel / smoothness figures for JopRBF (normalized compact-support Wendland-C2 RBF).
#
# Generates two PNGs under docs/JopRBF/images/:
#   1. jopRBF_basis_2d.png       a single normalized basis, its C2 vs C0 kernel, compact support
#   2. jopRBF_smoothness_2d.png  why C2: linear (C0) interpolation kinks between nodes; Wendland (C2) is smooth
#
# Run (global environment has JetPack dev'd + PyPlot):
#   julia --startup-file=no docs/JopRBF/JopRBF_demo.jl
#
ENV["MPLBACKEND"] = "Agg"
using JetPack, Jets, LinearAlgebra, Random, Printf, PyPlot

Random.seed!(20260709)
const OUT = joinpath(@__DIR__, "images")
isdir(OUT) || mkpath(OUT)

# scattered node cloud (denser shallow) used for the basis figure
function make_nodes(nz, nx)
    pts = Vector{Float64}[]
    for x in range(1, nx, length = 12)
        push!(pts, [1.0, x]); push!(pts, [Float64(nz), x])
    end
    for z in range(1, nz, length = 9)
        push!(pts, [z, 1.0]); push!(pts, [z, Float64(nx)])
    end
    for _ = 1:120
        push!(pts, [1 + (nz - 1) * rand()^0.65, 1 + (nx - 1) * rand()])
    end
    reduce(hcat, pts)
end

# normalized (Shepard) RBF evaluation on the full grid, for a chosen kernel
function eval_normalized(nodes, c, nz, nx, delta, kernel)
    d = zeros(nz, nx)
    M = size(nodes, 2)
    for iz = 1:nz, ix = 1:nx
        num = 0.0; den = 0.0
        for j = 1:M
            r = sqrt((iz - nodes[1, j])^2 + (ix - nodes[2, j])^2) / delta
            w = kernel(r)
            num += w * c[j]; den += w
        end
        d[iz, ix] = den > 0 ? num / den : 0.0
    end
    d
end

wend_c2(r) = r < 1 ? (1 - r)^4 * (4r + 1) : 0.0     # Wendland C2 kernel (what JopRBF uses)
linear0(r) = r < 1 ? (1 - r) : 0.0                  # linear "tent" kernel -> C0 (piecewise linear)

# ---------------------------------------------------------------------------
# Figure 1: a single normalized basis + C2 vs C0 kernel + compact support
# ---------------------------------------------------------------------------
nz, nx = 160, 220
delta = 26.0
nodes = make_nodes(nz, nx)
M = size(nodes, 2)
A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, nz, nx), nodes; delta = delta)

jc = argmin(vec(sum((nodes .- [80.0, 110.0]) .^ 2, dims = 1)))
e = zeros(domain(A)); e[jc] = 1.0
basis = A * e
zc, xc = nodes[1, jc], nodes[2, jc]

figure(figsize = (13, 4.2))
subplot(1, 3, 1)
im = imshow(basis, aspect = "auto", cmap = "viridis")
scatter([xc], [zc], s = 45, c = "red", marker = "x")
xlim(1, nx); ylim(nz, 1)
title("normalized basis  A e_j  (compact bump)"); xlabel("x"); ylabel("z"); colorbar(im, fraction = 0.046)

subplot(1, 3, 2)
rr = range(0, 1.3, length = 400)
plot(rr, wend_c2.(rr), lw = 2.4, label = "Wendland C2  (1-r)^4(4r+1)")
plot(rr, linear0.(rr), lw = 1.8, ls = "--", label = "linear C0  (1-r)")
axvline(1.0, color = "k", lw = 0.7, alpha = 0.5)
xlabel("r = |x - ξ| / δ"); ylabel("φ(r)"); title("kernel shape: C2 vs C0"); legend(); grid(alpha = 0.3)

subplot(1, 3, 3)
row = clamp(round(Int, zc), 1, nz)
plot(1:nx, basis[row, :], lw = 2.0)
axvline(xc - delta, color = "r", ls = ":", lw = 1); axvline(xc + delta, color = "r", ls = ":", lw = 1)
xlabel("x"); ylabel("basis value"); title("slice at z=$(row): zero beyond ±δ"); grid(alpha = 0.3)
tight_layout()
savefig(joinpath(OUT, "jopRBF_basis_2d.png"), dpi = 120); close()

# ---------------------------------------------------------------------------
# Figure 2: why C2. Same coarse nodes and coefficients, two interpolation
# kernels. Linear weights give a piecewise-linear (C0) field that kinks between
# nodes; the Wendland kernel (what JopRBF uses) is smooth (C2).
# ---------------------------------------------------------------------------
gz, gx = 130, 180
h = 26.0
snodes = let pts = Vector{Float64}[]
    zs = collect(range(1.0, gz, step = h)); (gz - zs[end] > 1) && push!(zs, Float64(gz))
    xs = collect(range(1.0, gx, step = h)); (gx - xs[end] > 1) && push!(xs, Float64(gx))
    for z in zs, x in xs
        jz = (1 < z < gz) ? 0.18h * (rand() - 0.5) : 0.0
        jx = (1 < x < gx) ? 0.18h * (rand() - 0.5) : 0.0
        push!(pts, [z + jz, x + jx])
    end
    reduce(hcat, pts)
end
sM = size(snodes, 2)
sdelta = 1.15h
crand = rand(sM)
dC0 = eval_normalized(snodes, crand, gz, gx, sdelta, linear0)
dC2 = eval_normalized(snodes, crand, gz, gx, sdelta, wend_c2)
cmin = min(minimum(dC0), minimum(dC2)); cmax = max(maximum(dC0), maximum(dC2))
row = 66

figure(figsize = (15, 4.4))
subplot(1, 3, 1)
im = imshow(dC0, aspect = "auto", cmap = "viridis", vmin = cmin, vmax = cmax)
scatter(snodes[2, :], snodes[1, :], s = 12, c = "white", edgecolors = "k", linewidths = 0.3)
axhline(row, color = "r", lw = 1.5)
xlim(1, gx); ylim(gz, 1)
title("linear kernel → C0 field\n(kinks/creases between nodes)"); xlabel("x"); ylabel("z"); colorbar(im, fraction = 0.046)

subplot(1, 3, 2)
im = imshow(dC2, aspect = "auto", cmap = "viridis", vmin = cmin, vmax = cmax)
scatter(snodes[2, :], snodes[1, :], s = 12, c = "white", edgecolors = "k", linewidths = 0.3)
axhline(row, color = "r", lw = 1.5)
xlim(1, gx); ylim(gz, 1)
title("Wendland C2 kernel → C2 field\n(smooth; what JopRBF uses)"); xlabel("x"); ylabel("z"); colorbar(im, fraction = 0.046)

subplot(1, 3, 3)
plot(1:gx, dC0[row, :], lw = 1.8, ls = "--", color = "tab:orange", label = "C0 linear (kinks)")
plot(1:gx, dC2[row, :], lw = 2.6, color = "tab:blue", label = "C2 Wendland (smooth)")
xlabel("x"); ylabel("field value along red line"); title("transect at z=$(row)")
legend(); grid(alpha = 0.3)
tight_layout()
savefig(joinpath(OUT, "jopRBF_smoothness_2d.png"), dpi = 120); close()

# scaling note
let x = rand(domain(A)), y = rand(range(A))
    A * x; A' * y
    t1 = @elapsed for _ = 1:5; A * x; end
    t2 = @elapsed for _ = 1:5; A' * y; end
    @printf("timing on %dx%d, M=%d, %d threads: forward %.1f ms, adjoint %.1f ms\n",
        nz, nx, M, Threads.nthreads(), 1000t1 / 5, 1000t2 / 5)
end
println("wrote jopRBF_basis_2d.png and jopRBF_smoothness_2d.png in ", OUT)
