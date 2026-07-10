#=
Figure generation for JopRBF in a water-bottom "freeze" workflow.

Builds a 101 x 201 (nz x nx) velocity model with a sloping water bottom, parameterizes only the sediment
below the water bottom with a depth-tapered scattered RBF node cloud, and freezes the water column. Produces
docs/JopRBF/images/jopRBF_waterbottom_2d.png.

Node placement is done by the hand-rolled `tapered_rbf_nodes` (rbf_tapered_nodes.jl, pure base Julia, NO
external mesher): brick-offset rows below the water bottom with spacing(z) = v_trend(z) / (freq * ppw(z)),
ppw linear in depth (dense shallow, sparse deep). The figure compares the parameterization at several
inversion frequencies (lower frequency = longer wavelength = coarser node spacing). Only PyPlot lives in
this directory's local project (docs/JopRBF/Project.toml), not in JetPack itself.

The reduced parameterization operator is  P = S_below ∘ A_RBF , where A_RBF maps node coefficients to the
fine grid and S_below (a JopDiagonal mask) restricts the update to below the water bottom. The full model is
    v(c) = v_water · (water mask)  +  P c
so the water column is frozen (never a function of c).

Run (uses the local project alongside this script):
  julia --startup-file=no --project=docs/JopRBF -t 8 docs/JopRBF/JopRBF_waterbottom_demo.jl
=#
ENV["MPLBACKEND"] = "Agg"
using JetPack, Jets, LinearAlgebra, Random, Printf, PyPlot
include(joinpath(@__DIR__, "rbf_tapered_nodes.jl"))

Random.seed!(20260709)
const OUT = joinpath(@__DIR__, "images")
isdir(OUT) || mkpath(OUT)

nz, nx = 101, 201
vwater = 1500.0
dz = dx = 10.0          # grid spacing [m]
ppw_top = 6.0           # points per wavelength at the water bottom (dense, shallow)
ppw_bot = 4.0           # points per wavelength at the model bottom (sparse, deep)
support = 1.6           # node support radius = support * local node spacing
node_ms = 22            # node marker size in the figure [points^2]
freqs = (6.0, 3.0, 1.5) # inversion frequencies to compare (lower freq -> coarser node spacing) [Hz]
flabel(f) = @sprintf("%g", f)   # frequency label without trailing zeros (6.0 -> "6", 1.5 -> "1.5")

# sloping water bottom: 12 samples deep on the left, 38 on the right
wbdepth(ix) = 12.0 + (38.0 - 12.0) * (ix - 1) / (nx - 1)

# "true" model: constant water above the water bottom, v(z) sediment plus two smooth anomalies below it
function true_model()
    v = fill(vwater, nz, nx)
    for ix = 1:nx, iz = 1:nz
        wb = wbdepth(ix)
        if iz > wb
            s = 1650 + 3.2 * (iz - wb)
            s += 180 * exp(-(((iz - 55) / 16)^2 + ((ix - 70) / 28)^2))
            s -= 140 * exp(-(((iz - 82) / 18)^2 + ((ix - 150) / 34)^2))
            v[iz, ix] = s
        end
    end
    v
end

vtrue = true_model()

# per-column water bottom (fine-grid index), below-WB mask, frozen water field, and a linear v(z) trend
wbcol = [clamp(round(Int, wbdepth(ix)), 1, nz) for ix = 1:nx]
below = [iz >= wbcol[ix] for iz = 1:nz, ix = 1:nx]
water_field = vwater .* .!below
a0, b0 = linear_vtrend(vtrue, wbcol)
vtrend = iz -> a0 + b0 * (iz - 1)

# build the RBF parameterization at a given frequency: tapered node cloud -> operator -> fit coefficients.
# precondition = false keeps the raw partition-of-unity kernel (best for a least-squares FIT / representation
# demo; the column-balancing precondition=true default is for reduced-parameterization inversion instead).
function build_rbf(freq)
    nodes, deltas = tapered_rbf_nodes(nz, nx; wb = wbcol, vtrend = vtrend, freq = freq, dz = dz,
        ppw_top = ppw_top, ppw_bot = ppw_bot, support = support)
    M = size(nodes, 2)
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, nz, nx), nodes; delta = deltas, precondition = false)
    P = JopDiagonal(Float64.(below)) ∘ A            # freeze the water column
    pou = A * ones(domain(A))
    c = convert(Matrix, P) \ vec(vtrue .* below)    # fit coefficients to the true sediment
    vmodel = reshape(P * c, nz, nx) .+ water_field
    relerr = (vmodel .- vtrue) ./ vtrue
    rms = sqrt(sum(abs2, relerr[below]) / count(below))
    @info "$(flabel(freq)) Hz: $M nodes ($(round(Int, count(below) / M))x reduction), min coverage $(round(minimum(pou[below]), digits=3)), RMS $(round(rms, digits=4))"
    (; nodes, vmodel, relerr, rms, pou, M)
end

res = [build_rbf(f) for f in freqs]

# ---------------------------------------------------------------------------
# figure: 4 rows x 2 cols
#   true model | coverage ;  6 Hz RBF | 6 Hz error ;  3 Hz | error ;  1.5 Hz | error
# ---------------------------------------------------------------------------
ext = [0.5, nx + 0.5, nz + 0.5, 0.5]
wbline = [wbdepth(ix) for ix = 1:nx]
vmin, vmax = 1500.0, maximum(vtrue)
emax = maximum(maximum(abs, r.relerr[below]) for r in res)   # common error color scale

function modelpanel(pos, img, ttl; nodes = nothing)
    subplot(4, 2, pos)
    im = imshow(img, extent = ext, aspect = "auto", cmap = "viridis", vmin = vmin, vmax = vmax)
    plot(1:nx, wbline, "w-", lw = 1.2)
    nodes === nothing || scatter(nodes[2, :], nodes[1, :], s = node_ms, c = "white", edgecolors = "k", linewidths = 0.3)
    xlim(1, nx); ylim(nz, 1); xlabel("x"); ylabel("z"); title(ttl); colorbar(im, fraction = 0.046)
end

function errpanel(pos, relerr, ttl)
    subplot(4, 2, pos)
    img = copy(relerr); img[.!below] .= NaN
    im = imshow(img, extent = ext, aspect = "auto", cmap = "seismic", vmin = -emax, vmax = emax)
    plot(1:nx, wbline, "k-", lw = 0.8)
    xlim(1, nx); ylim(nz, 1); xlabel("x"); ylabel("z"); title(ttl)
    cb = colorbar(im, fraction = 0.046); cb.set_label("Δv / v")
end

figure(figsize = (12, 16))

# row 1: true model | coverage (finest frequency)
modelpanel(1, vtrue, "true model")
subplot(4, 2, 2)
covimg = copy(res[1].pou); covimg[.!below] .= NaN
imc = imshow(covimg, extent = ext, aspect = "auto", cmap = "magma", vmin = 0.999, vmax = 1.001)
plot(1:nx, wbline, "c-", lw = 1.2); xlim(1, nx); ylim(nz, 1); xlabel("x"); ylabel("z")
title("coverage below WB  A·1 ≈ 1  ($(flabel(freqs[1])) Hz, water frozen white)"); colorbar(imc, fraction = 0.046)

# rows 2-4: one frequency per row
for (i, f) in enumerate(freqs)
    r = res[i]
    modelpanel(2i + 1, r.vmodel, "$(flabel(f)) Hz RBF  ($(r.M) nodes, ppw $(flabel(ppw_top))→$(flabel(ppw_bot)))"; nodes = r.nodes)
    errpanel(2i + 2, r.relerr, @sprintf("%s Hz error  (RMS %.4f)", flabel(f), r.rms))
end

tight_layout()
savefig(joinpath(OUT, "jopRBF_waterbottom_2d.png"), dpi = 120); close()
println("wrote ", joinpath(OUT, "jopRBF_waterbottom_2d.png"))
