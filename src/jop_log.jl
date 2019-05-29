"""
    F = JopLog(spc)

where `F` is the log operator with domain and range given by `spc::JetSpace`.
"""

JopLog(spc::JetSpace{T,N}) where {T,N} = JopNl(dom = spc, rng = spc, f! = JopLog_f!, df! = JopLog_df!, df′! = JopLog_df′!, upstate! = JopLog_upstate!, s = (invu = zeros(spc),))
export JopLog

function JopLog_f!(d::AbstractArray{T,N}, u::AbstractArray{T,N}; kwargs...) where {T,N}
    d .= log.(u)
end

# TODO: assert denominator > 0 
function JopLog_upstate!(u::AbstractArray, s::NamedTuple) 
    s.invu .= 1 ./ u
end

function JopLog_df!(δd::AbstractArray, δm::AbstractArray; invu, kwargs...) 
    δd .= invu .* δm
end

function JopLog_df′!(δm::AbstractArray, δd::AbstractArray; invu, kwargs...) 
    δm .= conj(invu) .* δd
end
