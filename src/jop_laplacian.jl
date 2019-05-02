"""
    A = JopLaplacian(sp)

Build a 2D (d^2/dx^2+d^2/dy^2) or 3D (d^2/dx^2+d^2/dy^2+d^2/dz^2) Laplacian
operator with three point centered finite difference stencil, with
`sp::JetSpace` and `A::JopLaplacian`.

For example,
```julia
A = JopLaplacian(JetSpace(Float64,128,256))
m = rand(domain(A))
d = A*m
```
"""
JopLaplacian(spc::JetSpace) = JopLn(dom = spc, rng = spc, df! = JopLaplacian_df!)

export JopLaplacian

function JopLaplacian_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kwargs...) where {T}
    nz,nx=size(m)
    for ix=1:nx, iz=1:nz
        d[iz,ix] = -4*m[iz,ix]
        if ix > 1
            d[iz,ix] += m[iz,ix-1]
        end
        if ix < nx
            d[iz,ix] += m[iz,ix+1]
        end

        if iz > 1
            d[iz,ix] += m[iz-1,ix]
        end
        if iz < nz
            d[iz,ix] += m[iz+1,ix]
        end
    end
    for ix=1:nx
        d[1,ix] += m[1,ix]
        d[nz,ix] += m[nz,ix]
    end
    for iz=1:nz
        d[iz,1] += m[iz,1]
        d[iz,nx] += m[iz,nx]
    end
    d
end

function JopLaplacian_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; kwargs...) where {T}
    nz,ny,nx=size(m)
    for ix=1:nx,iy=1:ny,iz=1:nz
        d[iz,iy,ix] = -6*m[iz,iy,ix]
        if ix > 1
            d[iz,iy,ix] += m[iz,iy,ix-1]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end
        if ix < nx
            d[iz,iy,ix] += m[iz,iy,ix+1]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end

        if iy > 1
            d[iz,iy,ix] += m[iz,iy-1,ix]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end
        if iy < ny
            d[iz,iy,ix] += m[iz,iy+1,ix]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end

        if iz > 1
            d[iz,iy,ix] += m[iz-1,iy,ix]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end
        if iz < nz
            d[iz,iy,ix] += m[iz+1,iy,ix]
        else
            d[iz,iy,ix] += m[iz,iy,ix]
        end
    end
    d
end
