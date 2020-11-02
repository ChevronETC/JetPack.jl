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
JopLaplacian(spc::JetSpace) = JopLn(dom = spc, rng = spc, df! = JopLaplacian_df!, df′! = JopLaplacian_df′!)

export JopLaplacian

function JopLaplacian_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kwargs...) where {T}
    nz,nx=size(m)
    d .= 0
    for ix=2:nx-1, iz=2:nz-1
        d[iz,ix] -= 4 * m[iz+0,ix+0]
        d[iz,ix] += m[iz-1,ix+0]
        d[iz,ix] += m[iz+1,ix+0]
        d[iz,ix] += m[iz+0,ix-1]
        d[iz,ix] += m[iz+0,ix+1]
    end
    d
end

function JopLaplacian_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; kwargs...) where {T}
    nz,nx=size(m)
    m .= 0
    for ix=2:nx-1, iz=2:nz-1
        m[iz+0,ix+0] -= 4 * d[iz,ix]
        m[iz-1,ix+0] += d[iz,ix]
        m[iz+1,ix+0] += d[iz,ix]
        m[iz+0,ix-1] += d[iz,ix]
        m[iz+0,ix+1] += d[iz,ix]
    end
    m
end

function JopLaplacian_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; kwargs...) where {T}
    nz,ny,nx=size(m)
    d .= 0
    for ix=2:nx-1,iy=2:ny-1,iz=2:nz-1
        d[iz,iy,ix] -= 6 * m[iz,iy,ix]
        d[iz,iy,ix] += m[iz-1,iy+0,ix+0]
        d[iz,iy,ix] += m[iz+1,iy+0,ix+0]
        d[iz,iy,ix] += m[iz+0,iy-1,ix+0]
        d[iz,iy,ix] += m[iz+0,iy+1,ix+0]
        d[iz,iy,ix] += m[iz+0,iy+0,ix-1]
        d[iz,iy,ix] += m[iz+0,iy+0,ix+1]
    end
    d
end

function JopLaplacian_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,3}; kwargs...) where {T}
    nz,ny,nx=size(m)
    m .= 0
    for ix=2:nx-1,iy=2:ny-1,iz=2:nz-1
        m[iz+0,iy+0,ix+0] -= 6 * d[iz,iy,ix]
        m[iz-1,iy+0,ix+0] += d[iz,iy,ix]
        m[iz+1,iy+0,ix+0] += d[iz,iy,ix]
        m[iz+0,iy-1,ix+0] += d[iz,iy,ix]
        m[iz+0,iy+1,ix+0] += d[iz,iy,ix]
        m[iz+0,iy+0,ix-1] += d[iz,iy,ix]
        m[iz+0,iy+0,ix+1] += d[iz,iy,ix]
    end
    m
end
