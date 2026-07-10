#=
rbf_tapered_nodes.jl -- generate a depth-tapered scattered node cloud for JopRBF, with NO external mesher.
Pure base Julia. Rows CONFORM to the water bottom: each row sits at a fixed offset from the (per-column) WB
and the offsets march downward by the local spacing, so the cloud is dense near the water bottom and coarser
at depth (no horizontal rows cross-cutting a dipping WB). Node spacing follows a points-per-wavelength rule,
spacing(z) = vtrend(z) / (freq * ppw(z)), with ppw LINEAR in depth from `ppw_top` at the water bottom to
`ppw_bot` at the model bottom. Rows are brick-offset (every other row shifted half a spacing) for near-
hexagonal packing; a few rows sit just above the WB (a capped "ghost" collar), a flat collar sits just below
the model bottom, and the lateral edges x=1/x=nx are pinned in every row, so all four boundaries get two-
sided support and the fit does not overshoot at any edge.
=#

"""
    nodes, deltas = tapered_rbf_nodes(nz, nx; wb, vtrend, freq, dz, ppw_top=5.0, ppw_bot=2.0, support=1.6, nghost=1, ghost_cap=4.0, nghost_bot=1)

Depth-tapered scattered RBF node cloud on an `nz x nx` grid, in fine-grid index coordinates, with rows that
CONFORM to (run parallel to) the water bottom. `wb` is a per-column `Vector` of fine-grid depth indices
(length `nx`) or a scalar; `vtrend(iz)` returns a reference velocity [m/s] at depth index `iz`; `dz` is the
grid spacing [m]; `support` sets each node's compact support radius as a multiple of its local spacing. Each
row sits at a fixed offset from the (per-column) water bottom (`z = wb[ix] + offset`); the offsets march
downward by the local spacing `vtrend(z) / (freq * ppw(z)) / dz`, so the cloud is dense at the water bottom
and coarser at depth (ppw LINEAR in depth from `ppw_top` at the WB to `ppw_bot` at the model bottom).
`nghost` rows sit just ABOVE the water bottom (each offset capped at `ghost_cap` fine samples, so they stay
near the WB even at low frequency) to give the near-WB fit two-sided support and prevent boundary overshoot.
`nghost_bot` FLAT rows sit just below the model bottom (mirror of the WB ghost) and the lateral edges
x=1/x=nx are pinned in every row, so all four boundaries get two-sided support. Rows are brick-offset (every
other row shifted half a spacing) for near-hexagonal packing, clipped to the grid. Returns `nodes` (2 x M,
row 1 = z, row 2 = x) and `deltas` (2 x M per-node isotropic support radii), both in fine-grid index units,
ready for `JopRBF(...; delta=deltas)`.
"""
function tapered_rbf_nodes(nz::Integer, nx::Integer; wb, vtrend, freq::Real, dz::Real,
                           ppw_top::Real = 5.0, ppw_bot::Real = 2.0, support::Real = 1.6,
                           nghost::Integer = 1, ghost_cap::Real = 4.0, nghost_bot::Integer = 1)
    wbcol = wb isa AbstractVector ? collect(Int, wb) : fill(Int(wb), nx)
    length(wbcol) == nx || error("tapered_rbf_nodes -- wb vector must have length nx = $nx")
    wbmin = minimum(wbcol)
    ppw(z)     = ppw_top + (ppw_bot - ppw_top) * clamp((z - (wbmin - 1)) / (nz - wbmin), 0.0, 1.0)
    spacing(z) = vtrend(z) / (freq * ppw(z)) / dz             # node spacing in FINE SAMPLES at depth index z

    # row offsets from the water bottom (fine samples): `nghost` capped rows ABOVE the WB (two-sided support),
    # the WB row (offset 0), then march DOWN by the local spacing to the model bottom.
    offs = Float64[]
    o = 0.0
    for _ = 1:nghost
        s = spacing(max(wbmin + o, 1.0)); s > 0 || error("tapered_rbf_nodes -- spacing must be positive; got $s")
        o -= min(s, ghost_cap); pushfirst!(offs, o)
    end
    push!(offs, 0.0)
    o = 0.0
    while wbmin + o <= nz
        s = spacing(wbmin + o); s > 0 || error("tapered_rbf_nodes -- spacing must be positive; got $s")
        o += s; wbmin + o <= nz && push!(offs, o)
    end

    # brick-offset lateral positions for a row of spacing s, always pinning the x=1 and x=nx boundaries so the
    # left/right edges get two-sided support (no one-sided lateral-edge overshoot).
    rowx(s, k) = begin
        xs = Float64[]; x = isodd(k) ? 1.0 : 1.0 + 0.5s
        while x <= nx; push!(xs, x); x += s; end
        (isempty(xs) || xs[1] > 1.0) && pushfirst!(xs, 1.0)
        xs[end] < nx && push!(xs, Float64(nx)); xs
    end

    # place brick-offset rows PARALLEL to the water bottom: node depth = wbcol[ix] + offset, clipped to grid.
    zp = Float64[]; xp = Float64[]; rp = Float64[]
    for (k, off) in enumerate(offs)
        s = spacing(wbmin + max(off, 0.0))                   # lateral spacing at this row's reference depth
        for x in rowx(s, k)
            ix = clamp(round(Int, x), 1, nx); z = wbcol[ix] + off
            1.0 <= z <= nz && (push!(zp, z); push!(xp, x); push!(rp, support * s))
        end
    end

    # bottom collar: `nghost_bot` FLAT rows just BELOW the model bottom (mirror of the WB ghost) so the deepest
    # cells get two-sided support instead of a one-sided edge (the model bottom is flat, so these rows are flat).
    s_bot = spacing(nz); zb = float(nz)
    for g = 1:nghost_bot
        zb += min(s_bot, ghost_cap)
        for x in rowx(s_bot, g + length(offs))
            push!(zp, zb); push!(xp, x); push!(rp, support * s_bot)
        end
    end

    nodes = permutedims(hcat(zp, xp))                        # (2, M): row 1 = z, row 2 = x
    nodes, permutedims(hcat(rp, rp))                         # per-node isotropic support radius (2, M)
end

"""
    a, b = linear_vtrend(vmodel, wbcol)

Least-squares linear fit `v(iz) = a + b*(iz-1)` [m/s] of the lateral-mean velocity below the water bottom,
a convenient `vtrend` for `tapered_rbf_nodes`. `wbcol` is the per-column water-bottom index.
"""
function linear_vtrend(vmodel::AbstractMatrix, wbcol::AbstractVector)
    nz, nx = size(vmodel)
    zc = Float64[]; vc = Float64[]
    for iz = 1:nz
        idx = [ix for ix = 1:nx if iz >= wbcol[ix]]
        isempty(idx) && continue
        push!(zc, iz - 1); push!(vc, sum(@view vmodel[iz, idx]) / length(idx))
    end
    a, b = hcat(ones(length(zc)), zc) \ vc
    a, b
end
