"""
    A = JopCubicSpline(dom, rng)

Cubic spline interpolation from coarse to fine grids.\\

'dom': operator domain. It corresponds the B-spline control points (nodes) on the coarse grid\\
'rng': operator range. It corresponds to the output array on the fine grid\\

The same constructor is used for 1D, 2D, 3D operators. The ND operator with N > 3 is simply a 3D operator applied to
every slice in the hyper dimensions. Currently, only the 4D is implemented.
"""
function JopCubicSpline(dom::JetSpace{T,N}, rng::JetSpace{T,N}) where {T,N}

    if N < 1
        error("Domain and range must be at least 1D")
    end

    n1 = size(rng)[1]
    n2 = (N > 1) ? size(rng)[2] : nothing
    n3 = (N > 2) ? size(rng)[3] : nothing
    nc1 = size(dom)[1]
    nc2 = (N > 1) ? size(dom)[2] : nothing
    nc3 = (N > 2) ? size(dom)[3] : nothing

    if nc1 > n1 || nc1 < 2
        error("In dimension 1, the domain size must be > 1 and <= range size")
    end

    if N > 1 && (nc2 > n2 || nc2 < 2)
        error("In dimension 2, the domain size must be > 1 and <= range size")
    end

    if N > 2 && (nc3 > n3 || nc3 < 2)
        error("In dimension 3, the domain size must be > 1 and <= range size")
    end

    if N > 3
        for i = 4:N
            if size(rng)[i] != size(dom)[i]
                error("Beyond dimension 3, domain and range must have the same size")
            end
        end
    end

    # compute the 1D kernel for each node using the basic B-spline function
    # also save the support minmax of each node on the dense grid
    d1 = convert(T, (n1 - 1) / (nc1 - 1))
    d2 = (N > 1) ? convert(T, (n2 - 1) / (nc2 - 1)) : nothing
    d3 = (N > 2) ? convert(T, (n3 - 1) / (nc3 - 1)) : nothing
        
    kernel1 = zeros(T, n1, nc1)
    kernel2 = (N > 1) ? zeros(T, n2, nc2) : nothing
    kernel3 = (N > 2) ? zeros(T, n3, nc3) : nothing
    iminmax1 = zeros(Int, nc1, 2)
    iminmax2 = (N > 1) ? zeros(Int, nc2, 2) : nothing
    iminmax3 = (N > 2) ? zeros(Int, nc3, 2) : nothing

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
    if N > 1
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
    if N > 2
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
        s = (; kernel1, kernel2, kernel3, iminmax1, iminmax2, iminmax3))
end

export JopCubicSpline

function JopCubicSpline_df!(d::AbstractArray{T,1}, m::AbstractArray{T,1}; kernel1, iminmax1, kwargs...) where {T}
    d .= 0
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            d[i] += kernel1[i,ic] * m[ic]
        end
    end
    d
end

function JopCubicSpline_df′!(m::AbstractArray{T,1}, d::AbstractArray{T,1}; kernel1, iminmax1, kwargs...) where {T}
    m .= 0
    nc1 = size(m,1)
    for ic = 1:nc1
        for i = iminmax1[ic,1]:iminmax1[ic,2]
            m[ic] += kernel1[i,ic] * d[i]
        end
    end
    m
end

function JopCubicSpline_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2, kwargs...) where {T}
    d .= 0
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
    d
end

function JopCubicSpline_df′!(m::AbstractArray{T,2}, d::AbstractArray{T,2}; kernel1, iminmax1, kernel2, iminmax2, kwargs...) where {T}
    m .= 0
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
    m
end

function JopCubicSpline_df!(d::AbstractArray{T,3}, m::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T}
    d .= 0
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    for ic3 = 1:nc3
        for k = iminmax3[ic3,1]:iminmax3[ic3,2]
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
    d
end

function JopCubicSpline_df′!(m::AbstractArray{T,3}, d::AbstractArray{T,3}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T}
    m .= 0
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    for ic3 = 1:nc3
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
    m
end

function JopCubicSpline_df!(d::AbstractArray{T,4}, m::AbstractArray{T,4}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T}
    d .= 0
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    nc = size(m,4)
    for c = 1:nc
        for ic3 = 1:nc3
            for k = iminmax3[ic3,1]:iminmax3[ic3,2]
                for ic2 = 1:nc2
                    for j = iminmax2[ic2,1]:iminmax2[ic2,2]
                        for ic1 = 1:nc1
                            for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                                d[i,j,k,c] += kernel1[i,ic1] * kernel2[j,ic2] * kernel3[k,ic3] * m[ic1,ic2,ic3,c]
                            end
                        end
                    end
                end
            end
        end
    end
    d
end

function JopCubicSpline_df′!(m::AbstractArray{T,4}, d::AbstractArray{T,4}; kernel1, iminmax1, kernel2, iminmax2, kernel3, iminmax3, kwargs...) where {T}
    m .= 0
    nc1 = size(m,1)
    nc2 = size(m,2)
    nc3 = size(m,3)
    nc = size(m,4)
    for c = 1:nc
        for ic3 = 1:nc3
            for k = iminmax3[ic3,1]:iminmax3[ic3,2]
                for ic2 = 1:nc2
                    for j = iminmax2[ic2,1]:iminmax2[ic2,2]
                        for ic1 = 1:nc1
                            for i = iminmax1[ic1,1]:iminmax1[ic1,2]
                                m[ic1,ic2,ic3,c] += kernel1[i,ic1] * kernel2[j,ic2] * kernel3[k,ic3] * d[i,j,k,c]
                            end
                        end
                    end
                end
            end
        end
    end
    m
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