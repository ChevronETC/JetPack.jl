"""
    JopMix(s::JetAbstractSpace{T,N}, nmix::NTuple{N,Int})

2D spatial mix, and is used to add smoothness to a domain.  JopMix is ported
from the CVX frequency domain FWI tools in SeisSpace.

# Notes:
* This appears to over-lap in funtionality with JopRoughness.

# Examples:
## 2D
```julia
nz, nx = 11,11
spc = JetSpace(Float32, nz, nx)
A = JopMix(spc, (5,5))
y = A * rand(dom)
x = A' * rand(rng)
```
"""
function JopMix(spc::JetAbstractSpace{T,N}, nmix::NTuple{N,Int}) where {T,N}
    for i in nmix
        i < 0 && error("mix lengths must be >= 0 -- $(nmix)")
    end
    JopLn(dom = spc, rng = spc, df! = JopMix_df!, df′! = JopMix_df′!, s= (nmix=nmix,))
end

export JopMix

function JopMix_df!(rngvec::AbstractArray{T,2}, domvec::AbstractArray{T,2}; nmix, kwargs...) where {T}
    rngvec .= 0
    nz,nx = size(domvec)
    lenz = div(nmix[1], 2)
    lenx = div(nmix[2], 2)
    wgt = (2 * lenz + 1) * (2 * lenx + 1)

    for kz = 1:nz
        jz1 = clamp(kz - lenz, 1, nz)
        jz2 = clamp(kz + lenz, 1, nz)
        for kx = 1:nx
            jx1 = clamp(kx - lenx, 1, nx)
            jx2 = clamp(kx + lenx, 1, nx)
            for jx = jx1:jx2
                for jz = jz1:jz2
                    rngvec[kz,kx] += domvec[jz,jx]
                end
            end
        end
    end
    rngvec ./= wgt
    rngvec
end

function JopMix_df′!(domvec::AbstractArray{T,2}, rngvec::AbstractArray{T,2}; nmix, kwargs...) where {T}
    domvec .= 0
    nz,nx = size(domvec)
    lenz = div(nmix[1], 2)
    lenx = div(nmix[2], 2)
    wgt = (2 * lenz + 1) * (2 * lenx + 1)

    for kz = 1:nz
        jz1 = clamp(kz - lenz, 1, nz)
        jz2 = clamp(kz + lenz, 1, nz)
        for kx = 1:nx
            jx1 = clamp(kx - lenx, 1, nx)
            jx2 = clamp(kx + lenx, 1, nx)
            for jx = jx1:jx2
                for jz = jz1:jz2
                    domvec[kz,kx] += rngvec[jz,jx]
                end
            end
        end
    end
    domvec ./= wgt
    domvec
end
