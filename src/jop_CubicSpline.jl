"""
    A = JopCubicSpline(dom, rng[, n_active_dimensions]; x1=nothing, x2=nothing, x3=nothing)

Cubic spline interpolation from coarse to fine grids.\\

'dom': operator domain. It corresponds to the B-spline control points (nodes) on the coarse grid over the active dimensions, augmented with passive dimensions.\\
'rng': operator range. It corresponds to the output array on the fine grid over the active dimensions, augmented with passive dimensions.\\
'n_active_dimensions': number of active dimensions. It can only be 1, 2, or 3 but <= than the number of dimensions of the domain and range.
'x1', 'x2', 'x3': optional 1D arrays of fractional node locations in [0,1], for active dimensions 1, 2, and 3 respectively.
If omitted, implicit regular spacing is used (`range(0,1,length=nc)`).

The same constructor is used for `N`-dimensional operators. The ND operator with `N > n_active_dimensions` is simply a `n_active_dimensions`-D operator applied to
every slice in the hyper dimensions.

# Examples:

```julia
A = JopCubicSpline(JetSpace(Float64,10), JetSpace(Float64,10))
m = ones(domain(A))
d = A*m ## input and output are 1D and have the same size
```

```julia
A = JopCubicSpline(JetSpace(Float64,5,2), JetSpace(Float64,10,2), 1)
m = ones(domain(A))
d = A*m ## input and output are 2D and have the same size along the second dimension; B-spline interpolation is applied along the first dimension only
```

```julia
A = JopCubicSpline(JetSpace(Float64,5,5,5), JetSpace(Float64,10,20,30))
B = JopCubicSpline(JetSpace(Float64,5,5,5), JetSpace(Float64,10,20,30), 3) ## same as A
m = ones(domain(A))
d = A*m ## input and output are 3D; B-spline interpolation is applied along all three dimensions
```

```julia
A = JopCubicSpline(JetSpace(Float64,5,5,5,2,3), JetSpace(Float64,10,20,30,2,3))
B = JopCubicSpline(JetSpace(Float64,5,5,5,2,3), JetSpace(Float64,10,20,30,2,3), 3) ## same as A
m = ones(domain(A))
d = A*m ## input and output are 5D and have the same size along the 4th and 5th dimensions; B-spline interpolation is applied along the first three dimensions
```

```julia
x1 = [0.0, 0.15, 0.33, 0.7, 1.0]
x2 = [0.0, 0.2, 0.5, 0.8, 1.0]
A = JopCubicSpline(JetSpace(Float64,5,5), JetSpace(Float64,20,30), 2; x1=x1, x2=x2)
```
"""
function JopCubicSpline(dom::JetSpace{T,N}, rng::JetSpace{T,N}, na::Int=0; x1::Union{Nothing, AbstractVector{<:Real}}=nothing, x2::Union{Nothing, AbstractVector{<:Real}}=nothing, x3::Union{Nothing, AbstractVector{<:Real}} =nothing) where {T,N}

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

    # if x1, x2, x3 are nothing, build regular kernels, otherwise build irregular kernels
    # dimension 1
    kernel1, iminmax1 = _build_cubic_spline_kernel(T, n1, nc1, x1)

    # dimension 2
    if na > 1
        kernel2, iminmax2 = _build_cubic_spline_kernel(T, n2, nc2, x2)
    else
        kernel2, iminmax2 = nothing, nothing
    end

    # dimension 3
    if na > 2
        kernel3, iminmax3 = _build_cubic_spline_kernel(T, n3, nc3, x3)
    else
        kernel3, iminmax3 = nothing, nothing
    end

    JopLn(; dom, rng, df! = JopCubicSpline_df!, df′! = JopCubicSpline_df′!, 
        s = (; kernel1, kernel2, kernel3, iminmax1, iminmax2, iminmax3, na))
end

export JopCubicSpline

function JopCubicSpline_df!(d::AbstractArray{T,N}, m::AbstractArray{T,N}; na, kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T,N}
    d .= 0

    # Loop over all dimensions beyond the active dimensions
    extra_dims = size(d)[na+1:N]
    indices = CartesianIndices(extra_dims)

    if na == 3
        for I in indices
            idx = I.I
            JopCubicSpline3D_df!(@view(d[:,:,:,idx...]), @view(m[:,:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3)
        end
        return d
    elseif na == 2
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSpline2D_df!(@view(d[:,:,idx...]), @view(m[:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2)
        end
        return d
    else
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSpline1D_df!(@view(d[:,idx...]), @view(m[:,idx...]); kernel1, iminmax1)
        end
        return d
    end
end

function JopCubicSpline_df′!(m::AbstractArray{T,N}, d::AbstractArray{T,N}; na, kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T,N}
    m .= 0

    # Loop over all dimensions beyond the active dimensions
    extra_dims = size(d)[na+1:N]
    indices = CartesianIndices(extra_dims)

    if na == 3
        for I in indices
            idx = I.I
            JopCubicSpline3D_df′!(@view(m[:,:,:,idx...]), @view(d[:,:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3)
        end
        return m
    elseif na == 2
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSpline2D_df′!(@view(m[:,:,idx...]), @view(d[:,:,idx...]); kernel1, iminmax1, kernel2, iminmax2)
        end
        return m
    else        
        Threads.@threads :static for I in indices
            idx = I.I
            JopCubicSpline1D_df′!(@view(m[:,idx...]), @view(d[:,idx...]); kernel1, iminmax1)
        end
        return m
    end
end


function JopCubicSpline1D_df!(d::AbstractArray{T,1}, m::AbstractArray{T,1}; kernel1, iminmax1) where {T}
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            d[i] += kernel1[i,ic] * m[ic]
        end
    end
end

function JopCubicSpline1D_df′!(m::AbstractArray{T,1}, d::AbstractArray{T,1}; kernel1, iminmax1) where {T}
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            m[ic] += kernel1[i,ic] * d[i]
        end
    end
    m
end

function JopCubicSpline2D_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2) where {T}
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

function JopCubicSpline2D_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2) where {T}
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

function JopCubicSpline3D_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3) where {T}
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

function JopCubicSpline3D_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3) where {T}
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

"""
Build cubic spline kernel for regularly or irregularly spaced nodes
"""
function _build_cubic_spline_kernel(T, n::Int, nc::Int, x::Union{Nothing, AbstractVector{<:Real}})
    if x === nothing
        xT = T.(collect(LinRange(0,1,nc)))
        dx = fill(T((n - 1) / (nc - 1)), nc)
    else
        (length(x) == nc) || error("x must have length equal to nc")
        minimum(x) >= 0.0 || error("x values must lie in [0,1]")
        maximum(x) <= 1.0 || error("x values must lie in [0,1]")
        allunique(x) || error("x values must not be duplicated")
        xT = T.(x)
        sort!(xT)
        xT[1] = T(0)
        xT[end] = T(1)
        # Treat near-uniform explicit coordinates as exactly regular to avoid
        # float32 quantization changing support bounds compared to x === nothing.
        dx_frac = diff(Float64.(xT))
        is_nearly_uniform = maximum(abs.(dx_frac .- dx_frac[1])) <= eps(T)
        if is_nearly_uniform
            xT = T.(collect(LinRange(0,1,nc)))
            dx = fill(T((n - 1) / (nc - 1)), nc)
        else
            dx = (n - 1) .* _node_local_spacing(xT)
        end
    end
    xmid = (n - 1) .* xT

    kernel = zeros(T, n, nc)
    iminmax = zeros(Int, nc, 2)

    @inbounds for ic = 1:nc
        imin = Int(floor(xmid[ic] - 2 * dx[ic]))
        imax = Int(floor(xmid[ic] + 2 * dx[ic]))
        imin = max(0, imin)
        imax = min(n, imax)
        iminmax[ic,1] = imin + 1
        iminmax[ic,2] = imax
        for i = imin:imax-1
            kernel[i+1,ic] = _cubic_spline(i, xmid=xmid[ic], dx=dx[ic])
        end
    end

    # add the missing contributions from out-of-bounds implicit nodes
    mid = - dx[1]
    imax = Int(floor(mid + 2 * dx[1]))
    imax = min(n, imax)
    @inbounds for i = 0:imax-1
        kernel[i+1,1] += _cubic_spline(i, xmid=mid, dx=dx[1])
    end

    mid = xmid[end] + dx[end]
    imin = Int(floor(mid - 2 * dx[end]))
    imin = max(0, imin)
    @inbounds for i = imin:n-1
        kernel[i+1,end] += _cubic_spline(i, xmid=mid, dx=dx[end])
    end

    kernel, iminmax
end

"""
Cubic spline basis function
"""
function _cubic_spline(x::Real; xmid::Real = 0, dx::Real = 1)
    xmin = xmid - 2 * dx
    xmax = xmid + 2 * dx
    t = convert(typeof(xmid),0)
    if x <= xmin || x >= xmax
        return t
    else
        t = (x - xmin) / dx
        if t <= 1
            return t^3 / 6
        elseif  t <= 2
            return (-3*t^3 + 12*t^2 - 12*t + 4) / 6
        elseif  t <= 3
            return (3*t^3 - 24*t^2 + 60*t - 44) / 6
        else
            return (-t^3 + 12*t^2 - 48*t + 64) / 6
        end
    end
end

function _node_local_spacing(x::Vector{<:Real})
    nc = length(x)
    dx = zeros(eltype(x), nc)
    dx[1] = x[2] - x[1]
    dx[end] = x[end] - x[end-1]
    @inbounds for i = 2:nc-1
        dx[i] = 0.5 * (x[i+1] - x[i-1])
    end
    dx
end