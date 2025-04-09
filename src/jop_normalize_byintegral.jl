"""
    F = JopNormalizeByIntegral(spc, c)

where `F` is the 'normalize by maximuam value of integral' operator F(t) = x(t) / int_0^T [ x(t) dt ]
Note the integration is along the 1st (fastest) dimension.
"""
function JopNormalizeByIntegral(spc::JetSpace{T}; α=1.0) where {T}
    tmp1 = zeros(spc)
    tmp2 = zeros(spc)
    A = JopLeakyIntegration(spc, α=α)
    JopNl(dom = spc, rng = spc, f! = JopNormalizeByIntegral_f!, df! = JopNormalizeByIntegral_df!, df′! = JopNormalizeByIntegral_df′!, s = (tmp1=tmp1, tmp2=tmp2, A=A))
end
export JopNormalizeByIntegral

function JopNormalizeByIntegral_f!(d::AbstractArray, m::AbstractArray; tmp1, tmp2, A)
    d .= m ./ sum(m, dims=1)
end

# For normalized integral method:
#   f   = I[m] / I_t[m]
#   df  = [ I[m]' I_t[m] - I[m] I_t[m]' ] / I_t[m]^2
#   df  = (I_t[m] I[dm] - I_t[dm] I[m]) / (I_t[m])^2
#   dfᵀ = (I_t[m] Iᵀ[δd] - I_tᵀ[mᵀ Iᵀ δd]) / (I_t[m])^2
#
# For normalization only:
#   f   = m / I_t[m]
#   df  = [ m' I_t[m] - m I_t[m]' ] / I_t[m]^2
#   df  = ( I_t[m] dm - I_t[dm] m) / (I_t[m])^2
#   dfᵀ = (I_t[m] δd - I_tᵀ[mᵀ δd]) / (I_t[m])^2

function JopNormalizeByIntegral_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ, tmp1, tmp2, A) where {T}
    δd .= (sum(mₒ,dims=1) .* δm .- sum(δm,dims=1) .* mₒ) ./ sum(mₒ,dims=1).^2
end

function JopNormalizeByIntegral_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ, tmp1, tmp2, A) where {T}
    δm .= conj.((sum(mₒ,dims=1) .* conj.(δd) .- sum(mₒ .* conj.(δd),dims=1)) ./ sum(mₒ,dims=1).^2)
end
