    using Printf

"""
    F = JopNormalize(spc, c)

where `F` is the 'normalize by maximuam value of integral' operator F(t) = x(t) / int_0^T [ x(t) dt ]
Note the integration is along the 1st (fastest) dimension.
"""
function JopNormalize(spc::JetSpace{T}; α=1.0) where {T}
    tmp1 = zeros(spc)
    tmp2 = zeros(spc)
    A = JopLeakyIntegration(spc, α=α)
    JopNl(dom = spc, rng = spc, f! = JopNormalize_f!, df! = JopNormalize_df!, df′! = JopNormalize_df′!, s = (tmp1=tmp1, tmp2=tmp2, A=A))
end
export JopNormalize

# d(t) = m(t) / [int_0^T m(t)]
function JopNormalize_f!(d::AbstractArray, m::AbstractArray; tmp1, tmp2, A)
    d .= 0
    n1 = size(d,1)
    n2 = div(length(d),n1)
    _d = reshape(d,n1,n2)
    _m = reshape(m,n1,n2)
    tmp1 .= A * m
    for k2 = 1:n2
        d[:,k2] .= _m[:,k2] ./ tmp1[end,k2]
    end
    d 
end

# d(t)  = m(t) / [int_0^T m(t)]
# d'(t) = m(t) (1 / [int_0^T m(t)])' + 1 / [int_0^T m(t)]

# (-I_t[dm] I[m] + I_t[m] I[dm]) / (I_t[m])^2
#δd .= (-sum(δm,dims=1) .* cumsum(mₒ,dims=1) .+ sum(mₒ,dims=1) .* cumsum(δm,dims=1)) ./ sum(mₒ,dims=1).^2
function JopNormalize_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ, tmp1, tmp2, A) where {T}
    n1 = size(δd,1)
    n2 = div(length(δd),n1)
    _δd = reshape(δd,n1,n2)
    _δm = reshape(δm,n1,n2)
    tmp1 .= A * mₒ
    tmp2 .= A * δm
    for k2 = 1:n2
        _δd[:,k2] .= (-tmp2[end,k2] .* tmp1[:,k2] .+ tmp1[end,k2] .* tmp2[:,k2]) ./ tmp1[end,k2]^2
    end
    δd 
end

# (-I_tᵀ[mᵀIᵀ] + I_t[m] Iᵀ) δd / (I_t[m])^2
#δm .= conj.(( -sum( mₒ .* (A'*conj.(δd)) ,dims=1) .+ ( sum(mₒ,dims=1) .* (A'*conj.(δd)) ) ) ./ sum(mₒ,dims=1).^2)
function JopNormalize_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ, tmp1, tmp2, A) where {T}
    n1 = size(δd,1)
    n2 = div(length(δd),n1)
    _δd = reshape(δd,n1,n2)
    _δm = reshape(δm,n1,n2)

    d1 = zeros(domain(A))
    d2 = zeros(domain(A))

    tmp_Atcd = A' * conj.(δd)
    tmp1 .= A * mₒ
    tmp2 .= A * (mₒ .* tmp_Atcd)
    for k2 = 1:n2
        _δm[:,k2] .= conj.(( -1 .* tmp2[end,k2] .+ (tmp1[end,k2] .* tmp_Atcd[:,k2]) ) ./ tmp1[end,k2].^2)

        d1[:,k2] .= tmp2[end,k2]
    end

    d2 .= sum( mₒ .* (A'*conj.(δd)) ,dims=1)

    for k2 = 1:n2
        @info @sprintf(" ")
        @info @sprintf("k2,n2,d1; %3d %3d %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f", k2, n2, d1[:,k2]...)
        @info @sprintf("k2,n2,d2; %3d %3d %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f", k2, n2, d2[:,k2]...)
    end

    #δm .= conj.(( -sum( mₒ .* (A'*conj.(δd)) ,dims=1) .+ ( sum(mₒ,dims=1) .* (A'*conj.(δd)) ) ) ./ sum(mₒ,dims=1).^2)

    δm 
end
