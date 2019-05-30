"""
    F = JopLog(spc)

where `F` is the log operator with domain and range given by `spc::JetSpace`.
"""
JopLog(spc::JetSpace) = JopNl(dom = spc, rng = spc, f! = JopLog_f!, df! = JopLog_df!, df′! = JopLog_df′!)
export JopLog

# TODO: assert denominator in linearization is > 0

JopLog_f!(d::AbstractArray, m::AbstractArray) = d .= log.(m)
JopLog_df!(δd::AbstractArray, δm::AbstractArray; mₒ) = δd .= δm ./ mₒ
JopLog_df′!(δm::AbstractArray, δd::AbstractArray; mₒ) = δm .= δd ./ conj(mₒ)
