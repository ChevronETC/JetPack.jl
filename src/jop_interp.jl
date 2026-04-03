"""
    JopInterp(dom, rng)

Performs linear interpolation from dom::JetSpace to rng::JetSpace. It is often
used to reduce dimensionality in FWI. JopInterp is ported from the CVX frequency
domain FWI tools in SeisSpace. Currently supports 2D (bi-linear) and 3D (tri-linear) 
interpolation. 

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

## 3D
```julia
dom = JetSpace(Float32, 11, 11, 11)
rng = JetSpace(Float32, 23, 23, 23)
A = JopInterp(dom, rng)
y = A * rand(dom)
x = A' * rand(rng)
```
"""
function JopInterp(dom::JetAbstractSpace, rng::JetAbstractSpace)
    ndim = length(size(dom))
    if ndim != length(size(rng))
        error("dimensionality of domain ($(ndim)) and range ($(length(size(rng)))) do not agree!")
    end

    if ndim < 2 || ndim > 3
        error("JopInterp supports 2D and 3D, got $(ndim)D")
    end

    dsize = size(dom)
    rsize = size(rng)

    for i in 1:ndim
        if rsize[i] < dsize[i]
            error("domain size ($(dsize[i])) exceeds range size ($(rsize[i])) in dimension $i!")
        end
        if dsize[i] < 2
            error("domain must have at least 2 points in dimension $i, got $(dsize[i])")
        end
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

    if isadj
        domvec .= 0
    else
        rngvec .= 0
    end
    
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

function JopInterp_apply_3d(domvec::AbstractArray, rngvec::AbstractArray, dom, rng, isadj::Bool)
    nzdom, nxdom, nydom = size(dom)
    nzrng, nxrng, nyrng = size(rng)
    zmin, zmax = 0.0, 1.0
    xmin, xmax = 0.0, 1.0
    ymin, ymax = 0.0, 1.0
    dzdom, dxdom, dydom = (zmax-zmin)/(nzdom-1), (xmax-xmin)/(nxdom-1), (ymax-ymin)/(nydom-1)
    dzrng, dxrng, dyrng = (zmax-zmin)/(nzrng-1), (xmax-xmin)/(nxrng-1), (ymax-ymin)/(nyrng-1)
    dxdydzdom = dxdom * dydom * dzdom
    ϵ = 1000 * eps(eltype(domvec))
    wmin, wmax = 0-ϵ, 1+ϵ

    if isadj
        domvec .= 0
    else
        rngvec .= 0
    end

    for kyrng = 1:nyrng
        for kxrng = 1:nxrng
            for kzrng = 1:nzrng
                zz = zmin + dzrng * (kzrng-1)
                xx = xmin + dxrng * (kxrng-1)
                yy = ymin + dyrng * (kyrng-1)

                kzdom = clamp(floor(Int, (zz - zmin) / dzdom) + 1, 1, nzdom-1)
                kxdom = clamp(floor(Int, (xx - xmin) / dxdom) + 1, 1, nxdom-1)
                kydom = clamp(floor(Int, (yy - ymin) / dydom) + 1, 1, nydom-1)

                z1 = zmin + dzdom * (kzdom - 1)
                z2 = zmin + dzdom * kzdom
                x1 = xmin + dxdom * (kxdom - 1)
                x2 = xmin + dxdom * kxdom
                y1 = ymin + dydom * (kydom - 1)
                y2 = ymin + dydom * kydom

                # trilinear weights (8 corners)
                w1 = (x2-xx) * (y2-yy) * (z2-zz) / dxdydzdom
                w2 = (xx-x1) * (y2-yy) * (z2-zz) / dxdydzdom
                w3 = (x2-xx) * (yy-y1) * (z2-zz) / dxdydzdom
                w4 = (xx-x1) * (yy-y1) * (z2-zz) / dxdydzdom
                w5 = (x2-xx) * (y2-yy) * (zz-z1) / dxdydzdom
                w6 = (xx-x1) * (y2-yy) * (zz-z1) / dxdydzdom
                w7 = (x2-xx) * (yy-y1) * (zz-z1) / dxdydzdom
                w8 = (xx-x1) * (yy-y1) * (zz-z1) / dxdydzdom

                wsum = w1+w2+w3+w4+w5+w6+w7+w8
                if abs(wsum - 1) > ϵ
                    error("sum of weights minus 1 ($(wsum-1)) not close to zero")
                end

                if w1<wmin||w1>wmax||w2<wmin||w2>wmax||w3<wmin||w3>wmax||w4<wmin||w4>wmax||
                   w5<wmin||w5>wmax||w6<wmin||w6>wmax||w7<wmin||w7>wmax||w8<wmin||w8>wmax
                    error("weights in interpolation exceed range [$(wmin),$(wmax)]")
                end

                iz0, iz1 = kzdom, kzdom+1
                ix0, ix1 = kxdom, kxdom+1
                iy0, iy1 = kydom, kydom+1

                if isadj
                    v = rngvec[kzrng, kxrng, kyrng]
                    domvec[iz0, ix0, iy0] += w1 * v
                    domvec[iz0, ix1, iy0] += w2 * v
                    domvec[iz0, ix0, iy1] += w3 * v
                    domvec[iz0, ix1, iy1] += w4 * v
                    domvec[iz1, ix0, iy0] += w5 * v
                    domvec[iz1, ix1, iy0] += w6 * v
                    domvec[iz1, ix0, iy1] += w7 * v
                    domvec[iz1, ix1, iy1] += w8 * v
                else
                    rngvec[kzrng, kxrng, kyrng] += w1 * domvec[iz0, ix0, iy0]
                    rngvec[kzrng, kxrng, kyrng] += w2 * domvec[iz0, ix1, iy0]
                    rngvec[kzrng, kxrng, kyrng] += w3 * domvec[iz0, ix0, iy1]
                    rngvec[kzrng, kxrng, kyrng] += w4 * domvec[iz0, ix1, iy1]
                    rngvec[kzrng, kxrng, kyrng] += w5 * domvec[iz1, ix0, iy0]
                    rngvec[kzrng, kxrng, kyrng] += w6 * domvec[iz1, ix1, iy0]
                    rngvec[kzrng, kxrng, kyrng] += w7 * domvec[iz1, ix0, iy1]
                    rngvec[kzrng, kxrng, kyrng] += w8 * domvec[iz1, ix1, iy1]
                end
            end
        end
    end
    nothing
end

JopInterp_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; dom, rng, kwargs...) where {T} = begin JopInterp_apply_2d(m, d, dom, rng, false); d end
JopInterp_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; dom, rng, kwargs...) where {T} = begin  JopInterp_apply_2d(m, d, dom, rng, true); m end

JopInterp_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; dom, rng, kwargs...) where {T} = begin JopInterp_apply_3d(m, d, dom, rng, false); d end
JopInterp_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,3}; dom, rng, kwargs...) where {T} = begin JopInterp_apply_3d(m, d, dom, rng, true); m end
