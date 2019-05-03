"""
    JopInterp(dom, rng)

Performs linear interpolation from dom::JetSpace to rng::JetSpace. It is often
used to reduce dimensionality in FWI.  JopInterp is ported from the CVX frequency
domain FWI tools in SeisSpace.

# Notes
* It is required that domain and range have the same dimensionality.
* It is required that domain and range both have more than 2 points per dimension
* It is required that the domain is coarser than the range.
* It is assumed that domain and range have the same boundary (e.g. in 2D: xmin,xmax,zmin,zmax).

# Examples:

## 2D
```julia
dom = JetSpace(Float32, 11, 11)
rng = JetSpace(Float32, 23, 23)
A = JopInterp(dom, rng)
y = A * rand(dom)
x = A' * rand(rng)
```
"""
function JopInterp(dom::JetAbstractSpace, rng::JetAbstractSpace)
    if length(size(dom)) != length(size(rng))
        error("dimensionality of domain ($(length(size(dom)))) and range ($(length(size(rng)))) do not agree!")
    end

    dsize = size(dom)
    rsize = size(rng)

    if rsize[1] < dsize[1] || rsize[2] < dsize[2]
        error("size of domain ($(size(dom))) smaller than range ($(size(rng))) !")
    end

    JopLn(dom=dom, rng=rng, df! = JopInterp_df!, df′! = JopInterp_df′!, s = (dom=dom, rng=rng))
end

export JopInterp

# coarse to dense
function JopInterp_apply_2d(domvec::AbstractArray, rngvec::AbstractArray, dom, rng, isadj::Bool)
    nzdom,nxdom = size(dom)
    nzrng,nxrng = size(rng)
    zmin,zmax = 0.0,1.0
    xmin,xmax = 0.0,1.0
    dzdom,dxdom = (zmax-zmin)/(nzdom-1),(xmax-xmin)/(nxdom-1)
    dzrng,dxrng = (zmax-zmin)/(nzrng-1),(xmax-xmin)/(nxrng-1)
    dxdzdom = dxdom * dzdom
    ϵ = 1000 * eps(eltype(domvec))
    wmin,wmax = 0-ϵ,1+ϵ

    for kxrng = 1:nxrng
        for kzrng = 1:nzrng
            xx = xmin + dxrng * (kxrng-1)
            zz = zmin + dzrng * (kzrng-1)

            kxdom = clamp(floor(Int, (xx - xmin) / dxdom) + 1, 1, nxdom-1)
            kzdom = clamp(floor(Int, (zz - zmin) / dzdom) + 1, 1, nzdom-1)

            x1 = xmin + dxdom * (kxdom - 1 + 0)
            x2 = xmin + dxdom * (kxdom - 1 + 1)
            z1 = zmin + dzdom * (kzdom - 1 + 0)
            z2 = zmin + dzdom * (kzdom - 1 + 1)

            w1 = (x2 - xx) * (z2 - zz) / dxdzdom
            w2 = (xx - x1) * (z2 - zz) / dxdzdom
            w3 = (x2 - xx) * (zz - z1) / dxdzdom
            w4 = (xx - x1) * (zz - z1) / dxdzdom

            # Sanity checks on weights
            if abs(w1 + w2 + w3 + w4 - 1) > ϵ
                error("sum of weights minus 1 ($(w1+w2+w3+w4-1)) not close to zero")
            end

            if w1 < wmin || w1 > wmax || w2 < wmin || w2 > wmax ||
                w3 < wmin || w3 > wmax || w4 < wmin || w4 > wmax
                error("weights in interpolation ($(w1),$(w2),$(w3),$(w1)) exceed range [$(wmin),$(wmax)]")
            end

            if isadj
                domvec[kzdom+0,kxdom+0] += w1 * rngvec[kzrng,kxrng]
                domvec[kzdom+0,kxdom+1] += w2 * rngvec[kzrng,kxrng]
                domvec[kzdom+1,kxdom+0] += w3 * rngvec[kzrng,kxrng]
                domvec[kzdom+1,kxdom+1] += w4 * rngvec[kzrng,kxrng]
            else
                rngvec[kzrng,kxrng] += w1 * domvec[kzdom+0,kxdom+0]
                rngvec[kzrng,kxrng] += w2 * domvec[kzdom+0,kxdom+1]
                rngvec[kzrng,kxrng] += w3 * domvec[kzdom+1,kxdom+0]
                rngvec[kzrng,kxrng] += w4 * domvec[kzdom+1,kxdom+1]
            end
        end
    end
    nothing
end

JopInterp_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; dom, rng, kwargs...) where {T} = begin JopInterp_apply_2d(m, d, dom, rng, false); d end
JopInterp_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; dom, rng, kwargs...) where {T} = begin  JopInterp_apply_2d(m, d, dom, rng, true); m end
