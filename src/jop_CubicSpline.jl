"""
    A = JopCubicSpline(dom, rng[, n_active_dimensions])

Cubic spline interpolation from coarse to fine grids.\\

'dom': operator domain. It corresponds to the B-spline control points (nodes) on the coarse grid over the active dimensions, augmented with passive dimensions.\\
'rng': operator range. It corresponds to the output array on the fine grid over the active dimensions, augmented with passive dimensions.\\
'n_active_dimensions': number of active dimensions. It can only be 1, 2, or 3 but <= than the number of dimensions of the domain and range.

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
"""
function JopCubicSpline(dom::JetSpace{T,N}, rng::JetSpace{T,N}, na::Int=0) where {T,N}

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

    # compute the 1D kernel for each node using the basic B-spline function
    # also save the support minmax of each node on the dense grid
    d1 = convert(T, (n1 - 1) / (nc1 - 1))
    d2 = (na > 1) ? convert(T, (n2 - 1) / (nc2 - 1)) : nothing
    d3 = (na > 2) ? convert(T, (n3 - 1) / (nc3 - 1)) : nothing
        
    kernel1 = zeros(T, n1, nc1)
    kernel2 = (na > 1) ? zeros(T, n2, nc2) : nothing
    kernel3 = (na > 2) ? zeros(T, n3, nc3) : nothing
    iminmax1 = zeros(Int, nc1, 2)
    iminmax2 = (na > 1) ? zeros(Int, nc2, 2) : nothing
    iminmax3 = (na > 2) ? zeros(Int, nc3, 2) : nothing

    # dimension 1
    for ic = 0:nc1-1
        mid = ic * d1
        imin = Int(floor(mid - 2 * d1))
        imax = Int(floor(mid + 2 * d1))
        imin = max(0, imin)
        imax = min(n1, imax)
        iminmax1[ic+1,1] = imin + 1
        iminmax1[ic+1,2] = imax
        for i = imin:imax-1
            kernel1[i+1,ic+1] = cubic_spline(i, xmid=mid, dx=d1)
        end
    end

    # add the missing contributions from out-of-bounds implicit nodes
    ic = 0
    mid = (ic - 1) * d1
    imax = Int(floor(mid + 2 * d1))
    imax = min(n1, imax)
    for i = 0:imax-1
        kernel1[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d1)
    end
    
    ic = nc1-1
    mid = (ic + 1) * d1
    imin = Int(floor(mid - 2 * d1))
    imin = max(0, imin)
    for i = imin:n1-1
        kernel1[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d1)
    end

    # dimension 2
    if na > 1
        for ic = 0:nc2-1
            mid = ic * d2
            imin = Int(floor(mid - 2 * d2))
            imax = Int(floor(mid + 2 * d2))
            imin = max(0, imin)
            imax = min(n2, imax)
            iminmax2[ic+1,1] = imin + 1
            iminmax2[ic+1,2] = imax
            for i = imin:imax-1
                kernel2[i+1,ic+1] = cubic_spline(i, xmid=mid, dx=d2)
            end
        end

        # add the missing contributions from out-of-bounds implicit nodes
        ic = 0
        mid = (ic - 1) * d2
        imax = Int(floor(mid + 2 * d2))
        imax = min(n2, imax)
        for i = 0:imax-1
            kernel2[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d2)
        end
        
        ic = nc2-1
        mid = (ic + 1) * d2
        imin = Int(floor(mid - 2 * d2))
        imin = max(0, imin)
        for i = imin:n2-1
            kernel2[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d2)
        end
    end
    
    # dimension 3
    if na > 2
        for ic = 0:nc3-1
            mid = ic * d3
            imin = Int(floor(mid - 2 * d3))
            imax = Int(floor(mid + 2 * d3))
            imin = max(0, imin)
            imax = min(n3, imax)
            iminmax3[ic+1,1] = imin + 1
            iminmax3[ic+1,2] = imax
            for i = imin:imax-1
                kernel3[i+1,ic+1] = cubic_spline(i, xmid=mid, dx=d3)
            end
        end

        # add the missing contributions from out-of-bounds implicit nodes
        ic = 0
        mid = (ic - 1) * d3
        imax = Int(floor(mid + 2 * d3))
        imax = min(n3, imax)
        for i = 0:imax-1
            kernel3[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d3)
        end
        
        ic = nc3-1
        mid = (ic + 1) * d3
        imin = Int(floor(mid - 2 * d3))
        imin = max(0, imin)
        for i = imin:n3-1
            kernel3[i+1,ic+1] += cubic_spline(i, xmid=mid, dx=d3)
        end
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
Uniform cubic spline basis function
"""
function cubic_spline(x::Real; xmid::Real = 0, dx::Real = 1)
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