# Figure generation for JopRBF in a water-bottom "freeze" workflow.
#
# Builds a 101 x 201 (nz x nx) velocity model with a sloping water bottom,
# parameterizes only the sediment below the water bottom with a scattered RBF node
# cloud, and freezes the water column. Produces
# docs/JopRBF/images/jopRBF_waterbottom_2d.png.
#
# Node placement is done by Gmsh (a JLL-backed mesh generator): the sediment below
# the water bottom is meshed with a size field spacing(z) = lambda(z) / ppw(z),
# where ppw(z) is linear in depth (denser shallow, sparser deep). The water bottom
# is embedded as a boundary curve, so mesh nodes land on it; the mesh vertices are
# used as the RBF centers. The figure compares the parameterization at several
# inversion frequencies (lower frequency = longer wavelength = coarser node
# spacing). Gmsh and PyPlot live in this directory's local project
# (docs/JopRBF/Project.toml), not in JetPack itself.
#
# The reduced parameterization operator is  P = S_below ∘ A_RBF , where A_RBF maps
# node coefficients to the fine grid and S_below (a JopDiagonal mask) restricts the
# update to below the water bottom. The full model is
#     v(c) = v_water · (water mask)  +  P c
# so the water column is frozen (never a function of c).
#
# Run (uses the local project alongside this script):
#   julia --startup-file=no --project=docs/JopRBF -t 8 docs/JopRBF/JopRBF_waterbottom_demo.jl
#
ENV["MPLBACKEND"] = "Agg"
using JetPack, Jets, LinearAlgebra, Random, Printf, PyPlot
using Gmsh: gmsh

Random.seed!(20260709)
const OUT = joinpath(@__DIR__, "images")
isdir(OUT) || mkpath(OUT)

nz, nx = 101, 201
vwater = 1500.0
dz = dx = 10.0          # grid spacing [m]
ppw_zmin = 6.0          # points per wavelength at zmin (the water bottom)
ppw_zmax = 4.0          # points per wavelength at zmax (model bottom); == ppw_zmin => resolution-matched (spacing tracks wavelength). Lower it to deliberately coarsen with depth.
support = 1.5           # node support radius = support * local spacing
max_deep_ratio = 1.35   # cap the coarsest node spacing at this * the shallow spacing (fills deep/bottom holes; Inf = no cap)
mesh_smoothing = 100    # Gmsh Laplacian smoothing passes ("go hard" on iterative refinement)
node_ms = 25            # node marker size in the figure [points^2]
freqs = (6.0, 3.0, 1.5) # inversion frequencies to compare (lower freq -> coarser node spacing) [Hz]
flabel(f) = @sprintf("%g", f)   # frequency label without trailing zeros (6.0 -> "6", 1.5 -> "1.5")

# sloping water bottom: 12 samples deep on the left, 38 on the right
wbdepth(ix) = 12.0 + (38.0 - 12.0) * (ix - 1) / (nx - 1)

# "true" model: constant water above the water bottom, v(z) sediment plus two
# smooth anomalies below it
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

# Gmsh meshing of the sediment below the water bottom with a depth-dependent size
# field spacing(z) = v_avg(z) / (freq * ppw(z)), ppw linear from ppw_zmin (WB) to
# ppw_zmax (model bottom). Lower freq => longer wavelength => coarser node spacing.
# The WB is embedded so mesh nodes land on it. Returns nodes (2 x M) and per-node
# support radii (2 x M).
function gmsh_nodes(freq)
    minwb = minimum(wbdepth(ix) for ix = 1:nx)
    rfun(z) = (1650.0 + 0.4 * max(0.0, z - minwb) * dz) /
              (freq * (ppw_zmin + (ppw_zmax - ppw_zmin) * clamp((z - minwb) / (nz - minwb), 0.0, 1.0))) / dz

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("wb")

    # domain = trapezoid below the water bottom, corners in (x, z) = (lateral, depth)
    p1 = gmsh.model.geo.addPoint(1.0, wbdepth(1), 0.0)          # WB, left
    p2 = gmsh.model.geo.addPoint(Float64(nx), wbdepth(nx), 0.0) # WB, right
    p3 = gmsh.model.geo.addPoint(Float64(nx), Float64(nz), 0.0) # bottom right
    p4 = gmsh.model.geo.addPoint(1.0, Float64(nz), 0.0)         # bottom left
    lwb = gmsh.model.geo.addLine(p1, p2)                        # water bottom (embedded => nodes on it)
    lr = gmsh.model.geo.addLine(p2, p3)
    lb = gmsh.model.geo.addLine(p3, p4)
    ll = gmsh.model.geo.addLine(p4, p1)
    cl = gmsh.model.geo.addCurveLoop([lwb, lr, lb, ll])
    gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()

    # size field: spacing(depth y) = vavg(y) / (fmean * ppw(y)) / dz, ppw linear in y
    frac = "max(0,min(1,(y-$(minwb))/$(nz - minwb)))"
    ppwe = "($(ppw_zmin)+($(ppw_zmax)-$(ppw_zmin))*$(frac))"
    vavg = "(1650+0.4*max(0,y-$(minwb))*$(dz))"
    sizeexpr = "$(vavg)/($(freq)*$(ppwe))/$(dz)"
    fid = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(fid, "F", sizeexpr)
    gmsh.model.mesh.field.setAsBackgroundMesh(fid)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    # cap the coarsest spacing (relative to the shallow spacing) so the deep/bottom
    # does not develop gaps; scales with frequency so all panels are treated alike
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_deep_ratio * 1650.0 / (freq * ppw_zmin) / dz)
    gmsh.option.setNumber("Mesh.Algorithm", 6)                 # Frontal-Delaunay (uniform, well-shaped triangles)
    gmsh.option.setNumber("Mesh.Smoothing", mesh_smoothing)    # iterative Laplacian smoothing of interior nodes
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(2)

    _, coord, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()

    c = reshape(coord, 3, :)                                    # rows: x, y(=z), z(=0)
    xs = c[1, :]; zs = c[2, :]
    nodes = permutedims(hcat(zs, xs))                          # (2, M): row1 = z, row2 = x
    rads = Float64[support * rfun(z) for z in zs]
    deltas = permutedims(hcat(rads, rads))
    nodes, deltas
end

vtrue = true_model()

# below-water-bottom mask and the frozen water field (shared across frequencies)
below = [iz > wbdepth(ix) for iz = 1:nz, ix = 1:nx]
water_field = vwater .* .!below

# build the RBF parameterization at a given frequency: mesh -> operator -> fit
function build_rbf(freq)
    nodes, deltas = gmsh_nodes(freq)
    M = size(nodes, 2)
    A = JopRBF(JetSpace(Float64, M), JetSpace(Float64, nz, nx), nodes; delta = deltas)
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
#   true model | coverage ;  6 Hz RBF | 6 Hz error ;  3 Hz | error ;  1 Hz | error
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
    modelpanel(2i + 1, r.vmodel, "$(flabel(f)) Hz RBF  ($(r.M) nodes)"; nodes = r.nodes)
    errpanel(2i + 2, r.relerr, @sprintf("%s Hz error  (RMS %.4f)", flabel(f), r.rms))
end

tight_layout()
savefig(joinpath(OUT, "jopRBF_waterbottom_2d.png"), dpi = 120); close()
println("wrote ", joinpath(OUT, "jopRBF_waterbottom_2d.png"))
