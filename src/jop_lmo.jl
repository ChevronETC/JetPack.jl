"""
    A = JopLMO(spc [;v=2.0, dx=1.0, dt=1.0])

2D linear moveout operator for domain and range `sp::JetSpace`.  The forward operator
maps from flat events to dipping events using the parameters `v,dx,dt`, and:

    (z_moveout) = (z + (dx/dt)/v * x)

# Notes

* It should be trivial to add a 3D implementation for this operator.
"""
JopLMO(spc::JetSpace{T}, v, dx, dt) where {T} = JopLN(dom = spc, rng = spc, df! = JopLMO_df!, df′! = JopLMO_df′!, s = (v=T((dx/dt)/v),))
JopLMO(spc::JetSpace{T}, v=2.0) where {T} = JopLn(dom = spc, rng = spc, df! = JopLMO_df!, df′! = JopLMO_df′!, s = (v=T(v),))

export JopLMO

function JopLMO_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; v, kwargs...) where {T}
    nz, nx = size(m)
    d .= 0
    for ix = 1:nx, iz = 1:nz
        iz_mo = iz + v*ix
        iz_mo_lb = floor(Int, iz_mo)
        iz_mo_ub = ceil(Int, iz_mo)
        if 1 <= iz_mo_lb <= nz && 1 <= iz_mo_ub <= nz
            a = iz_mo_ub == iz_mo_lb ? 0.0 : (iz_mo - iz_mo_lb) / (iz_mo_ub - iz_mo_lb)
            d[iz,ix] = (1 - a) * m[iz_mo_lb,ix] + a * m[iz_mo_ub,ix]
        end
    end
    d
end

function JopLMO_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; v, kwargs...) where {T}
    nz, nx = size(d)
    fill!(m, 0.0)
    for ix = 1:nx, iz = 1:nz
        iz_mo = iz + v*ix
        iz_mo_lb = floor(Int, iz_mo)
        iz_mo_ub = ceil(Int, iz_mo)
        if 1 <= iz_mo_lb <= nz && 1 <= iz_mo_ub <= nz
            a = iz_mo_ub == iz_mo_lb ? 0.0 : (iz_mo - iz_mo_lb) / (iz_mo_ub - iz_mo_lb)
            m[iz_mo_lb,ix] += (1 - a) * d[iz,ix]
            m[iz_mo_ub,ix] += a * d[iz,ix]
        end
    end
    m
end
