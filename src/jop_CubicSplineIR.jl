"""
    A = JopCubicSplineIR(dom, rng[, n_active_dimensions]; x1=nothing, x2=nothing, x3=nothing)

Cubic spline interpolation from coarse to fine grids.\\

'dom': operator domain. It corresponds to the B-spline control points (nodes) on the coarse grid over the active dimensions, augmented with passive dimensions.\\
'rng': operator range. It corresponds to the output array on the fine grid over the active dimensions, augmented with passive dimensions.\\
'n_active_dimensions': number of active dimensions. It can only be 1, 2, or 3 but <= than the number of dimensions of the domain and range.
'x1', 'x2', 'x3': optional 1D arrays of fractional node locations in [0,1], strictly ascending, for active dimensions 1, 2, and 3 respectively.
If omitted, regular spacing is used (`range(0,1,length=nc)`).

The same constructor is used for `N`-dimensional operators. The ND operator with `N > n_active_dimensions` is simply a `n_active_dimensions`-D operator applied to
every slice in the hyper dimensions.

# Examples:

```julia
A = JopCubicSplineIR(JetSpace(Float64,10), JetSpace(Float64,10))
m = ones(domain(A))
d = A*m ## input and output are 1D and have the same size
```

```julia
A = JopCubicSplineIR(JetSpace(Float64,5,2), JetSpace(Float64,10,2), 1)
m = ones(domain(A))
d = A*m ## input and output are 2D and have the same size along the second dimension; B-spline interpolation is applied along the first dimension only
```

```julia
A = JopCubicSplineIR(JetSpace(Float64,5,5,5), JetSpace(Float64,10,20,30))
B = JopCubicSplineIR(JetSpace(Float64,5,5,5), JetSpace(Float64,10,20,30), 3) ## same as A
m = ones(domain(A))
d = A*m ## input and output are 3D; B-spline interpolation is applied along all three dimensions
```

```julia
A = JopCubicSplineIR(JetSpace(Float64,5,5,5,2,3), JetSpace(Float64,10,20,30,2,3))
B = JopCubicSplineIR(JetSpace(Float64,5,5,5,2,3), JetSpace(Float64,10,20,30,2,3), 3) ## same as A
m = ones(domain(A))
d = A*m ## input and output are 5D and have the same size along the 4th and 5th dimensions; B-spline interpolation is applied along the first three dimensions
```

```julia
x1 = [0.0, 0.15, 0.33, 0.7, 1.0]
x2 = [0.0, 0.2, 0.5, 0.8, 1.0]
A = JopCubicSplineIR(JetSpace(Float64,5,5), JetSpace(Float64,20,30), 2; x1=x1, x2=x2)
```
"""
function _regular_fractional_locations_ir(nc::Int)
    collect(LinRange(0.0, 1.0, nc))
end

function _validate_fractional_locations_ir(x::AbstractVector{<:Real}, nc::Int, dim::Int)
    length(x) == nc || error("x$(dim) must have length equal to domain size in active dimension $(dim)")
    nc >= 2 || error("In active dimension $(dim), domain size must be > 1")

    x1 = Float64(x[1])
    xn = Float64(x[end])
    (x1 >= 0.0 && xn <= 1.0) || error("x$(dim) values must lie in [0,1]")

    @inbounds for i = 2:nc
        Float64(x[i]) > Float64(x[i-1]) || error("x$(dim) must be strictly ascending")
    end
    nothing
end

function _node_local_spacing_ir(xf::Vector{Float64})
    nc = length(xf)
    dx = zeros(Float64, nc)
    dx[1] = xf[2] - xf[1]
    dx[end] = xf[end] - xf[end-1]
    @inbounds for i = 2:nc-1
        dx[i] = 0.5 * (xf[i+1] - xf[i-1])
    end
    dx
end

function _build_cubic_spline_kernel_ir(T, n::Int, nc::Int, xf_in::AbstractVector{<:Real})
    xf = Float64.(xf_in)
    xmid = (n - 1) .* xf
    dxloc = (n - 1) .* _node_local_spacing_ir(xf)

    kernel = zeros(T, n, nc)
    iminmax = zeros(Int, nc, 2)

    @inbounds for ic = 1:nc
        mid = xmid[ic]
        dx = dxloc[ic]
        imin = Int(floor(mid - 2 * dx))
        imax = Int(floor(mid + 2 * dx))
        imin = max(0, imin)
        imax = min(n, imax)
        iminmax[ic,1] = imin + 1
        iminmax[ic,2] = imax
        for i = imin:imax-1
            kernel[i+1,ic] = cubic_spline(i, xmid=mid, dx=dx)
        end
    end

    mid = xmid[1] - dxloc[1]
    dx = dxloc[1]
    imax = Int(floor(mid + 2 * dx))
    imax = min(n, imax)
    @inbounds for i = 0:imax-1
        kernel[i+1,1] += cubic_spline(i, xmid=mid, dx=dx)
    end

    mid = xmid[end] + dxloc[end]
    dx = dxloc[end]
    imin = Int(floor(mid - 2 * dx))
    imin = max(0, imin)
    @inbounds for i = imin:n-1
        kernel[i+1,end] += cubic_spline(i, xmid=mid, dx=dx)
    end

    kernel, iminmax
end

function _build_uniform_cubic_spline_kernel_ir(T, n::Int, nc::Int)
    d = convert(T, (n - 1) / (nc - 1))
    kernel = zeros(T, n, nc)
    iminmax = zeros(Int, nc, 2)

    for ic = 0:nc-1
        mid = ic * d
        imin = Int(floor(mid - 2 * d))
        imax = Int(floor(mid + 2 * d))
        imin = max(0, imin)
        imax = min(n, imax)
        iminmax[ic+1,1] = imin + 1
        iminmax[ic+1,2] = imax
        for i = imin:imax-1
            kernel[i+1,ic+1] = cubic_spline(i, xmid=mid, dx=d)
        end
    end

    ic = 0
    mid = (ic - 1) * d
    imax = Int(floor(mid + 2 * d))
    imax = min(n, imax)
    for i = 0:imax-1
        kernel[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d)
    end

    ic = nc - 1
    mid = (ic + 1) * d
    imin = Int(floor(mid - 2 * d))
    imin = max(0, imin)
    for i = imin:n-1
        kernel[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d)
    end

    kernel, iminmax
end

function JopCubicSplineIR(dom::JetSpace{T,N}, rng::JetSpace{T,N}, na::Int=0; x1=nothing, x2=nothing, x3=nothing) where {T,N}

    if N < 1
        error("Domain and range must be at least 1D")
    end

    (na <= 0) && (na = min(3,N))

    for i = na+1:N
        if size(dom)[i] != size(rng)[i]
            error("In passive dimensions, domain and range must have the same size")
        end
    end

    # active dimensions
    n1 = size(rng)[1]
    n2 = (na > 1) ? size(rng)[2] : nothing
    n3 = (na > 2) ? size(rng)[3] : nothing
    nc1 = size(dom)[1]
    nc2 = (na > 1) ? size(dom)[2] : nothing
    nc3 = (na > 2) ? size(dom)[3] : nothing

    if nc1 > n1 || nc1 < 2
        error("In active dimension 1, the domain size must be > 1 and <= range size")
    end

    if na > 1 && (nc2 > n2 || nc2 < 2)
        error("In active dimension 2, the domain size must be > 1 and <= range size")
    end

    if na > 2 && (nc3 > n3 || nc3 < 2)
        error("In active dimension 3, the domain size must be > 1 and <= range size")
    end

    if na == 1 && x2 !== nothing
        error("x2 was provided but n_active_dimensions = 1")
    end
    if na <= 2 && x3 !== nothing
        error("x3 was provided but n_active_dimensions <= 2")
    end

    kernel1, iminmax1 = if x1 === nothing
        _build_uniform_cubic_spline_kernel_ir(T, n1, nc1)
    else
        _validate_fractional_locations_ir(x1, nc1, 1)
        _build_cubic_spline_kernel_ir(T, n1, nc1, x1)
    end

    kernel2, iminmax2 = if na > 1
        if x2 === nothing
            _build_uniform_cubic_spline_kernel_ir(T, n2, nc2)
        else
            _validate_fractional_locations_ir(x2, nc2, 2)
            _build_cubic_spline_kernel_ir(T, n2, nc2, x2)
        end
    else
        nothing, nothing
    end

    kernel3, iminmax3 = if na > 2
        if x3 === nothing
            _build_uniform_cubic_spline_kernel_ir(T, n3, nc3)
        else
            _validate_fractional_locations_ir(x3, nc3, 3)
            _build_cubic_spline_kernel_ir(T, n3, nc3, x3)
        end
    else
        nothing, nothing
    end

    JopLn(; dom, rng, df! = JopCubicSplineIR_df!, df′! = JopCubicSplineIR_df′!, 
        s = (; kernel1, kernel2, kernel3, iminmax1, iminmax2, iminmax3, na))
end

export JopCubicSplineIR

function JopCubicSplineIR_df!(d::AbstractArray{T,N}, m::AbstractArray{T,N}; na, kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T,N}
    d .= 0

    # Loop over all dimensions beyond the active dimensions
    extra_dims = size(d)[na+1:N]
    indices = CartesianIndices(extra_dims)

    if na == 3
        for I in indices
            idx = I.I
            JopCubicSplineIR3D_df!(@view(d[:,:,:,idx...]), @view(m[:,:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3)
        end
        return d
    elseif na == 2
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSplineIR2D_df!(@view(d[:,:,idx...]), @view(m[:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2)
        end
        return d
    else
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSplineIR1D_df!(@view(d[:,idx...]), @view(m[:,idx...]); kernel1, iminmax1)
        end
        return d
    end
end

function JopCubicSplineIR_df′!(m::AbstractArray{T,N}, d::AbstractArray{T,N}; na, kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T,N}
    m .= 0

    # Loop over all dimensions beyond the active dimensions
    extra_dims = size(d)[na+1:N]
    indices = CartesianIndices(extra_dims)

    if na == 3
        for I in indices
            idx = I.I
            JopCubicSplineIR3D_df′!(@view(m[:,:,:,idx...]), @view(d[:,:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3)
        end
        return m
    elseif na == 2
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSplineIR2D_df′!(@view(m[:,:,idx...]), @view(d[:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2)
        end
        return m
    else        
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSplineIR1D_df′!(@view(m[:,idx...]), @view(d[:,idx...]); kernel1, iminmax1)
        end
        return m
    end
end


function JopCubicSplineIR1D_df!(d::AbstractArray{T,1}, m::AbstractArray{T,1}; kernel1, iminmax1) where {T}
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            d[i] += kernel1[i,ic] * m[ic]
        end
    end
end

function JopCubicSplineIR1D_df′!(m::AbstractArray{T,1}, d::AbstractArray{T,1}; kernel1, iminmax1) where {T}
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            m[ic] += kernel1[i,ic] * d[i]
        end
    end
    m
end

function JopCubicSplineIR2D_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2) where {T}
    nc1 = size(m,1)
    nc2 = size(m,2)
    for ic2 = 1:nc2
        for j = iminmax2[ic2,1]:iminmax2[ic2,2]
            for ic1 = 1:nc1
                for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                    d[i,j] += kernel1[i,ic1] * kernel2[j,ic2] * m[ic1,ic2]
                end
            end
        end
    end
end

function JopCubicSplineIR2D_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2) where {T}
    nc1 = size(m,1)
    nc2 = size(m,2)
    for ic2 = 1:nc2
        for j = iminmax2[ic2,1]:iminmax2[ic2,2]
            for ic1 = 1:nc1
                for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                    m[ic1,ic2] += kernel1[i,ic1] * kernel2[j,ic2] * d[i,j]
                end
            end
        end
    end
end

function JopCubicSplineIR3D_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3) where {T}
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    for ic3 = 1:nc3
        Threads.@threads :static for k = iminmax3[ic3,1]:iminmax3[ic3,2]
            for ic2 = 1:nc2
                for j = iminmax2[ic2,1]:iminmax2[ic2,2]
                    for ic1 = 1:nc1
                        for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                            d[i,j,k] += kernel1[i,ic1] * kernel2[j,ic2] * kernel3[k,ic3] * m[ic1,ic2,ic3]
                        end
                    end
                end
            end
        end
    end
end

function JopCubicSplineIR3D_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3) where {T}
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    Threads.@threads :static for ic3 = 1:nc3
        for k = iminmax3[ic3,1]:iminmax3[ic3,2]
            for ic2 = 1:nc2
                for j = iminmax2[ic2,1]:iminmax2[ic2,2]
                    for ic1 = 1:nc1
                        for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                            m[ic1,ic2,ic3] += kernel1[i,ic1] * kernel2[j,ic2] * kernel3[k,ic3] * d[i,j,k]
                        end
                    end
                end
            end
        end
    end
end