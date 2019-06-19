"""
    F = JopNim(spc, c)

where `F` is the 'normalized integral method' operator F = int_0^t [ x(t) dt ] / int_0_T [ x(t) dt ]
"""
JopNim(spc::JetSpace{T}) where {T} = JopNl(dom = spc, rng = spc, f! = JopNim_f!, df! = JopNim_df!, df′! = JopNim_df′!)
export JopNim

function JopNim_f!(d::AbstractArray, m::AbstractArray)
    d .= cumsum(m,dims=1) ./ sum(m,dims=1)
end

# (-I_t[dm] I[m] + I_t[m] I[dm]) / (I_t[m])^2
function JopNim_df!(δd::AbstractArray, δm::AbstractArray; mₒ) 
    δd .= (-sum(δm,dims=1) .* cumsum(mₒ,dims=1).+ sum(mₒ,dims=1) .* cumsum(δm,dims=1)) ./ sum(mₒ,dims=1).^2
end

# (-I_tᵀ[mᵀIᵀ] + I_t[m] Iᵀ) δd / (I_t[m])^2
function JopNim_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ) where {T}
    A = JopLeakyIntegration(JetSpace(T,size(δm)), α=1.0)
    δm .= conj.(( -sum( mₒ .* (A'*conj.(δd)) ,dims=1) .+ ( sum(mₒ,dims=1) .* (A'*conj.(δd)) ) ) ./ sum(mₒ,dims=1).^2)
end
