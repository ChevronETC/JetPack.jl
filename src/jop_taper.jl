"""
    A = JopTaper(spc, dims, frac[, frac_end; mode=:normal, taper=(:cosine,:cosine)])

The linear operator `A` tapers the edges of all or some subset of the dimensions of an array that
belongs to `spc::JetSpace`.  The dimensions that are tapered are given by `dims`, an `Int` tuple.
`A` will taper the beginning and/or end edges of each specified array dimension. The size of the
beginning and end tapers are determined as a fraction of the length of that array dimension, and
this fraction is set using the tuple `frac`, and (optionally) `frac_end`.

If `frac_end::NTuple{M,Float64}` is not set, then the ith entry of `frac::NTuple{M,Float64}`
determines the length of the center portion of the ith dimension that is not tapered.  If
`frac_end::NTuple{M,NTuple{2,Float64}}` is set, then the length of the begining taper of the ith array
dimension is determiend by `frac[i]`, and the length of the end taper is determined by `frac_end[i]`.

The optional named argument `mode` can be used if the taper is applied to a dimension that
corresponds to FFT ordering where the edges are assumed to be at the center and left ends of the dimension.

The optional named argument `taper` is a tuple specifying what type of taper to use at each end.  Available
tapers are :cosine and :heviside.

# Examples:

## taper for a space containing 1D arrays
```julia
A = JopTaper(JetSpace(Float64,10), (1,), (.75,))
m = ones(domain(A))
d = A*m
```

## taper for a space containing 1D arrays, and only taper at the end
```julia
A = JopTaper(JetSpace(Float64,10), (1,), (0.0,), (0.25,))
m = ones(domain(A))
d = A*m
```

## taper for a space containing 2D arrays, where both dimensions are tapered
```julia
A = JopTaper(JetSpace(Float64,10,11), (1,2), (0.75,0.5))
m = ones(domain(A))
d = A*m
```

## taper for a space containing 3D arrays, where two of the three dimensions are tapered
```julia
A = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.75,0.5))
m = ones(domain(A))
d = A*m
```
In the next example, array dimension `1` is tapered at the end, and array dimension
`2` is tapered at the beginning and end.
```julia
A = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.0,0.25), (0.25,0.25))
m = ones(domain(A))
d = A*m
```
"""
function JopTaper(spc::JetAbstractSpace, dim::NTuple{M,Int}, prcbeg::NTuple{M,Real}, prcend::NTuple{M,Real}; mode=:normal, taper=(:cosine,:cosine)) where {M}
    # specify taper via middle part
    @assert length(prcbeg) == length(prcend) == length(dim)
    prc = [[Float64(prcbeg[i]),Float64(prcend[i])] for i=1:length(prcbeg)]
    mode = isa(mode,Symbol) ? ntuple(i->mode, length(dim)) : mode
    JopLn(dom = spc, rng = spc, df! = JopTaper_df!, s = (dim=dim, prc=prc, mode=mode, taper=taper))
end

function JopTaper(spc::JetAbstractSpace, dim::NTuple{M,Int}, prc::NTuple{M,Real}; mode=:normal, taper=(:cosine,:cosine)) where {M}
    # specify taper via middle part
    @assert length(dim) == length(prc)
    prc_lftrht = ntuple(idim->0.5*(1.0 - prc[idim]), length(dim))
    prc = [[Float64(prc_lftrht[idim]),Float64(prc_lftrht[idim])] for idim=1:length(dim)]
    mode = isa(mode,Symbol) ? ntuple(i->mode, length(dim)) : mode
    JopLn(dom = spc, rng = spc, df! = JopTaper_df!, s = (dim=dim, prc=prc, mode=mode, taper=taper))
end

export JopTaper

function JopTaper_df!(d::AbstractArray, m::AbstractArray; dim, prc, mode, taper, kwargs...)
    d .= m
    for idim = 1:length(dim)
        tpr = buildtaper(dim[idim], prc[idim], mode[idim], size(m, dim[idim]), taper)
        for idx in CartesianIndices(size(d))
            d[idx] *= tpr[idx[dim[idim]]]
        end
    end
    d
end

function buildtaper(dim, prc, mode, n, taper)
    tpr = ones(n)

    nlft = round(Int, prc[1] * n)
    nrht = round(Int, prc[2] * n)

    for i=1:nlft
        if taper[1] == :cosine
            x = (i - 1.0) / (nlft - 1.0)
            tpr[i] = 0.5 * (1.0 + cos(pi * (1.0 + x)))
        else # :heviside
            tpr[i] = 0.0
        end
    end

    for i=1:nrht
        if taper[2] == :cosine
            x = (i - 1.0) / (nrht - 1.0)
            tpr[n+1-i] = 0.5 * (1.0 + cos(pi * (1.0 + x)))
        else # :heviside
            tpr[n+1-i] = 0.0
        end
    end

    tpr = mod == :fft ? fftshift(tpr) : tpr
end
