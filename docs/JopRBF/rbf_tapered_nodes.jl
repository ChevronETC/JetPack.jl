#=
rbf_tapered_nodes.jl -- generate a depth-tapered scattered node cloud for JopRBF, with NO external mesher.
Pure base Julia. The node cloud is dense near the surface / water bottom and coarser at depth:
node spacing follows a points-per-wavelength rule, spacing(z) = vtrend(z) / (freq * ppw(z)), with ppw
LINEAR in depth from `ppw_top` at the water bottom to `ppw_bot` at the model bottom. Rows are brick-offset
(every other row shifted half a spacing) for near-hexagonal packing, clipped to at/below a (per-column or
flat) water bottom, plus one water-bottom-conforming row so the partition-of-unity coverage reaches the WB.
=#

"""
    nodes, deltas = tapered_rbf_nodes(nz, nx; wb, vtrend, freq, dz, ppw_top=5.0, ppw_bot=2.0, support=1.6)

Depth-tapered scattered RBF node cloud on an `nz x nx` grid, in fine-grid index coordinates. `wb` is a
per-column `Vector` of fine-grid depth indices (length `nx`) or a scalar; `vtrend(iz)` returns a reference
velocity [m/s] at depth index `iz`; `dz` is the grid spacing [m]; `support` sets each node's compact support
radius as a multiple of its local spacing. Returns `nodes` (2 x M, row 1 = z, row 2 = x) and `deltas`
(2 x M per-node isotropic support radii), both in fine-grid index units, ready for `JopRBF(...; delta=deltas)`.
"""
function tapered_rbf_nodes(nz::Integer, nx::Integer; wb, vtrend, freq::Real, dz::Real,
                           ppw_top::Real = 5.0, ppw_bot::Real = 2.0, support::Real = 1.6)
    wbcol = wb isa AbstractVector ? collect(Int, wb) : fill(Int(wb), nx)
    length(wbcol) == nx || error("tapered_rbf_nodes -- wb vector must have length nx = $nx")
    wbmin = minimum(wbcol)
    ppw(z)     = ppw_top + (ppw_bot - ppw_top) * clamp((z - (wbmin - 1)) / (nz - wbmin), 0.0, 1.0)
    spacing(z) = vtrend(z) / (freq * ppw(z)) / dz              # node spacing in FINE SAMPLES at depth index z

    zp = Float64[]; xp = Float64[]
    z = float(wbmin); parity = 0                               # brick-offset rows, WB -> model bottom
    while z <= nz
        s = spacing(z); s > 0 || error("tapered_rbf_nodes -- spacing must be positive; got $s at z=$z")
        x = iseven(parity) ? 1.0 : 1.0 + 0.5s                 # shift every other row by half a spacing
        while x <= nx
            ix = clamp(round(Int, x), 1, nx)
            z >= wbcol[ix] && (push!(zp, z); push!(xp, x))    # keep only at/below the per-column water bottom
            x += s
        end
        z += s; parity += 1
    end
    x = 1.0                                                    # one WB-conforming row for coverage up to the WB
    while x <= nx
        ix = clamp(round(Int, x), 1, nx); zwb = float(wbcol[ix])
        push!(zp, zwb); push!(xp, x); x += spacing(zwb)
    end

    nodes = permutedims(hcat(zp, xp))                          # (2, M): row 1 = z, row 2 = x
    rads  = Float64[support * spacing(z) for z in zp]
    nodes, permutedims(hcat(rads, rads))                       # per-node isotropic support radius (2, M)
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
