"""
    A = JopHighpass(sp)

Build a 2D Highpass operator, which first apply a center weighted three point
smoothing operator [0.25, 0.5, 0.25] to each dimension with a given iterations
number (nit[1:2]) to get the low frequency component of the data, then subtract
it from the original data to return the high pass component.  It works with
`sp::JetSpace`, `nit::Array{Integer, 1}` and `A::JopHighpass`.  It can specify
different number of iterations (nit) in x and z direction. For a seismic image,
more iterations in z direction is often applied to obtain more vertical
resolution.

For example,
```julia
A = JopHighpass(JetSpace(Float64,128,256), nit=(3,1)])
m = rand(domain(A))
d = A*m
```
"""
JopHighpass(spc::JetSpace; nit=(1,1)) = JopLn(dom = spc, rng = spc, df! = JopHighpass_df!, s = (nit=nit,))

export JopHighpass

function JopHighpass_df!(d::AbstractArray{T,2}, m::AbstractArray{T,2}; nit, kwarg...) where {T}
    if all(nit .== 0)
        copyto!(d,m)
        return d
    end

    nz,nx = size(m)

    temp1 = zeros(T, nz+2, nx)
    temp2 = zeros(T, nz, nx)
    padding_z!(temp1, m)
    for it = 1:nit[1]
        for ix = 1:nx
            for iz = 1:nz
                temp2[iz, ix] = 0.25*(temp1[iz,ix] + temp1[iz+2,ix]) + 0.5 * temp1[iz+1,ix]
            end
        end
        padding_z!(temp1, temp2)
    end
    temp3 = zeros(T,nz,nx+2)
    temp4 = zeros(T,nz,nx)
    padding_x!(temp3, temp1[2:nz+1,:])
    for it = 1:nit[2]
        for ix = 1:nx
            for iz = 1:nz
                temp4[iz, ix] = 0.25*(temp3[iz,ix] + temp3[iz,ix+2]) + 0.5 * temp3[iz,ix+1]
            end
        end
        padding_x!(temp3, temp4)
    end
    for ix = 1:nx
        for iz = 1:nz
            d[iz,ix] = m[iz,ix] - temp3[iz,ix+1]
        end
    end
    d
end

function padding_x!(dout::Array{T,2}, din::Array{T,2}) where T<:AbstractFloat
    dout[:,2:end-1] = din
    dout[:,1] = dout[:,2]
    dout[:,end] = dout[:,end-1]
end

function padding_z!(dout::Array{T,2}, din::Array{T,2}) where T<:AbstractFloat
    dout[2:end-1,:] = din
    dout[1,:] = dout[2,:]
    dout[end,:] = dout[end-1,:]
end
