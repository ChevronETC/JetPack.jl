"""
    A = JopPad(dom, pad..., [extend=false, accumulate=false])

where `dom::JetSpace` is the domain of `A`, and `pad::UnitRange...` determines the range of `A`.
If `extend=false`, then the padded region is set to zero.  If `extend=true`, then the padded
region is set from the boundary of the domain.  The `accumulate=true` option is specific to
the `Ginsu` operation in JetPackWave, and should not be used unless you really know what you
are doing.

# Examples:

## 1D
```julia
A = JopPad(JetSpace(Float64,2), -1:3)
m = [1.0, 2.0]
d = A*m # d = [0.0, 0.0, 1.0, 2.0, 0.0]
A = JopPad(JetSpace(Float64,2), -1:3, extend=true)
d = A*m # d = [1.0, 1.0, 1.0, 2.0, 2.0]
```

## 2D
```julia
A = JopPad(JetSpace(Float64,2,2), -1:3, 1:3)
m = [11. 12. ; 21. 22.]
d = A*m # d = [0. 0. 11. 12. 0. ; 0. 0. 21. 22. 0. ; 0. 0. 0. 0. 0.]
```

# Notes:
* This operator may also be used for truncation
* There may be overlap between the functionality of `JopPad` and `JopRestriction`, and it may be worth
thinking of how to consolidate them.
"""
function JopPad(dom::JetSpace{T}, pad::UnitRange...; extend=false, accumulate=false) where {T}
    ndims(dom) == length(pad) || error("'pad' options and the number of dimensions in the domain must be consistent")
    JopLn(dom = dom, rng = JetSpace(T, map(length, pad)), df! = JopPad_df!, df′! = JopPad_df′!,
        s = (pad=pad, extend=extend, accumulate=accumulate))
end
export JopPad

function JopPad_df!(d::AbstractArray, m::AbstractArray; pad, extend, kwargs...)
    indices_rng = CartesianIndices(size(d))
    indices_dom = CartesianIndices(size(m))
    mpadₒ = modeloffset(pad)
    for didx in indices_rng
        midx = shiftidx(mpadₒ, didx)
        if midx ∈ indices_dom
            d[didx] = m[midx]
        elseif extend == false
            d[didx] = 0.0
        elseif extend == true
            d[didx] = m[modelidx_nn(mpadₒ, didx, size(m))]
        end
    end
    d
end

function JopPad_df′!(m::AbstractArray{T,N}, d::AbstractArray{T,N}; accumulate, pad, extend, kwargs...) where {T,N}
    if accumulate
        JopPad_df′!(Val(extend), m, d, (x,y)->x+y, pad)
    else
        m .= zero(T)
        JopPad_df′!(Val(extend), m, d, (x,y)->y, pad)
    end
    m
end

function JopPad_df′!(extend::Val{true}, m::AbstractArray{T,N}, d::AbstractArray{T,N}, _::Function, pad) where {T,N}
    indices_rng = CartesianIndices(size(d))
    indices_dom = CartesianIndices(size(m))
    mpadₒ = modeloffset(pad)
    for didx in indices_rng
        midx = shiftidx(mpadₒ, didx)
        if midx ∈ indices_dom
            m[midx] += d[didx]
        else
            m[modelidx_nn(mpadₒ, didx, size(m))] += d[didx]
        end
    end
    nothing
end

function JopPad_df′!(extend::Val{false}, m::AbstractArray{T,N}, d::AbstractArray{T,N}, f::Function, pad) where {T,N}
    dpadₒ = dataoffset(pad)
    for midx in intersection_domain(pad, size(m))
        didx = shiftidx(dpadₒ, midx)
        m[midx] = f(m[midx], d[didx])
    end
    nothing
end

function intersection_domain(pad::NTuple{N}, n::NTuple{N}) where {N}
    function _intersection_domain(i, pad, n)
        j = clamp(pad[i][1], 1, n[i])
        k = clamp(pad[i][end], 1, n[i])
        j:k
    end
    CartesianIndices(ntuple(i->_intersection_domain(i,pad,n), Val(N)))
end

modeloffset(pad::NTuple{N,UnitRange}) where {N} = ntuple(i->pad[i][1]-1, Val(N))
dataoffset(pad::NTuple{N,UnitRange}) where {N} = ntuple(i->1-pad[i][1], Val(N))
shiftidx(padₒ::NTuple{N,Int}, idx::CartesianIndex{N}) where {N} = CartesianIndex{N}(ntuple(dim->idx[dim]+padₒ[dim], N))
function modelidx_nn(padₒ::NTuple{N,Int}, idx::CartesianIndex{N}, n::NTuple{N,Int}) where N
    midx_try = shiftidx(padₒ, idx)
    CartesianIndex{N}(ntuple(dim->clamp(midx_try[dim], 1, n[dim]), Val(N)))
end
