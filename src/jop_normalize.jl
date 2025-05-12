"""
    F = JopNormalize(spc(T), ϵ=eps(T), mode=:trace)

where `F` is the self-adjoint normalization operator F(t) = x(t) / || x(t) ||. 
Can be run in 'normalize by norm of trace' (mode=:trace), or 'normalize by norm of shot'
(mode=:shot) modes.  
Note the normalization is along the 1st dimension.

Copied from the implementation in Guillaume Barnier's SEP thesis
https://stacks.stanford.edu/file/druid:td173jf2299/Guillaume_Barnier_PhD_thesis-augmented.pdf
"""
function JopNormalize(spc::JetSpace{T}; ϵ::T = eps(T), mode::Symbol=:trace) where {T}
    JopNl(dom = spc, rng = spc, f! = JopNormalize_f!, df! = JopNormalize_df!, df′! = JopNormalize_df′!, s = (ϵ=ϵ, mode=mode, ))
end
export JopNormalize

function JopNormalize_f!(d::AbstractArray, m::AbstractArray; ϵ, mode)
    n1 = size(d,1) # Time index 
    n2 = size(d,2) # Trace index 
    n3 = div(length(d),n2*n1) # Component index

    d .= m 
    @views for k3 ∈ 1:n3
        (mode === :shot) && (norm = sqrt(dot(m[:,:,k3],m[:,:,k3])))
        for k2 ∈ 1:n2
            (mode === :trace) && (norm = sqrt(dot(m[:,k2,k3],m[:,k2,k3])))
            d[:,k2,k3] ./= (norm + ϵ)
        end
    end
    d
end

function JopNormalize_df!(δd::AbstractArray{T}, δm::AbstractArray{T}; mₒ, ϵ, mode) where {T}
    n1 = size(δd,1)
    n2 = size(δd,2)
    n3 = div(length(δd),n1*n2)
    
    _mₒ = reshape(mₒ,n1,n2,n3)
    _δd = reshape(δd,n1,n2,n3)
    _δm = reshape(δm,n1,n2,n3)
    δd .= 0

    @views for k3 ∈ 1:n3
        (mode === :shot) && (_corr = dot(_mₒ[:,:,k3], _δm[:,:,k3]))
        (mode === :shot) && (_norm = dot(_mₒ[:,:,k3], _mₒ[:,:,k3]))
        for k2 ∈ 1:n2
            (mode === :trace) && (_corr = dot(_mₒ[:,k2,k3], _δm[:,k2,k3]))
            (mode === :trace) && (_norm = dot(_mₒ[:,k2,k3], _mₒ[:,k2,k3]))
            _δd[:,k2,k3] .= _δm[:,k2,k3] ./ (sqrt(_norm) + ϵ) .- _mₒ[:,k2,k3] .* _corr ./ ((sqrt(_norm) + ϵ)^2 * sqrt(_norm + ϵ))
        end
    end
    δd
end

function JopNormalize_df′!(δm::AbstractArray{T}, δd::AbstractArray{T}; mₒ, ϵ, mode) where {T}
    JopNormalize_df!(δm, δd; mₒ, ϵ, mode)
end
