"""
    A = JopReghost(spc, zt, zm, dt, dx)

Build a 2D receiver-side reghosting and redatuming operator in the FK domain using
simple phase-shift operations (e.g. Posthumus, 1993).  `spc::JetSpaceSymmetric` is
the domain and range of `A`. `zt` is target depth where the receivers will be
re-dataumed to, and `zm` is the measurement depth where the recording is made.
`dt` is the time sampling interval and `dx` is the receiver sampling interval.

Since this operator is built in the **FK** domain, `spc` has symmetry in the
frequency dimension.  Hence `spc::JotSpaceSymmetric`.  We provide a convenience
method:

    sp = symspace(JopReghost, T, nw, nk)

where `T` is either `Float32` or `Float64`, `nw` are the number of frequency samples,
and `nk` are the number of wavenumber samples.

# Example:
```julia
P = JopPad(JetSpace(Float64,nt,nx), 1:2nt, 1:nx)
F = JopFft(range(P))
G = JopReghost(range(F), 20.0, 20.0, .004, 10.0)
A = P'*F'*G*F*P
data_reghost = A * data_noghost
```
Note that in the above example,
```julia
range(F) == symspace(JopReghost_df!, Float64, size(range(P))...)
```
"""
function JopReghost(spc::JetSSpace, zt, zm, dt, dx)
    nw = size(spc,1)
    nk = size(spc,2)
    nkpos = div(nk,2)+1
    nkneg = nk - nkpos
    wn = pi/dt
    dw = wn/nw
    kn = pi/dx
    dk = 2*kn/nk
    w = collect(range(0.0, stop=dw*(nw-1), length=nw))
    kpos = range(0.0, stop=dk*(nkpos-1), length=nkpos)
    kneg = range(-dk*nkneg, stop=-dk, length=nkneg)
    k = [kpos;kneg]

    @assert Base.length(w) == nw
    @assert Base.length(k) == nk
    @assert Base.length(kpos) == nkpos
    @assert Base.length(kneg) == nkneg

    JopLn(dom = spc, rng = spc, df! = JopReghost_df!, df′! = JopReghost_df′!, s = (zt=zt, zm=zm, dt=dt, dx=dx, k=k, w=w))
end

export JopReghost

JopReghost_df!(d::AbstractArray{Complex{T},2}, m::AbstractArray{Complex{T},2}; zt, zm, dt, dx, k, w, kwargs...) where {T} = reghost(d, m, -one(T), zt, zm, dt, dx, k, w)
JopReghost_df′!(m::AbstractArray{Complex{T},2}, d::AbstractArray{Complex{T},2}; zt, zm, dt, dx, k, w, kwargs...) where {T} = reghost(m, d, one(T), zt, zm, dt, dx, k, w)

function reghost(d::AbstractArray{Complex{T}}, m, sgn, zt, zm, dt, dx, k, w) where {T}
    nw, nk = size(d)
    @assert size(m) == (nw,nk)

    dzs = zt - zm
    dzg = zt + zm
    for iw = 1:nw
        wc = w[iw]^2 / T(1500)^2
        for ik = 1:nk
            kz2 = wc - k[ik]^2
            if kz2 >= 0.0
                kz = sgn*sqrt(kz2)
                d[iw,ik] = (exp(im*kz*dzs) - exp(im*kz*dzg)) * m[iw,ik]  # todo -- obliquity factor
            else
                d[iw,ik] = zero(T)
            end
        end
    end
    d
end

function JopReghost_symspace_map(I, n, _n)
    if I[1] > _n[1]
        return CartesianIndex((2 + n[1] - I[1], I[2]))
    else
        return CartesianIndex(I)
    end
end

function Jets.symspace(JopReghost, ::Type{T}, nw, nk) where {T}
    n = (nw,nk)
    _n = (div(nw,2)+1, nk)
    JetSSpace(Complex{T}, n, _n, I->JopReghost_symspace_map(I, n, _n))
end
