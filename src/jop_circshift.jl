"""
    A = JopCircShift(R::JetSpace, shifts)

`d=A*m` is equivalent to `circshift(d, m, shifts)`, and where `m` and `d` are in the domain/range of `A`.
"""
function JopCircShift(spc, shifts)
    JopLn(df! = JopCircShift_df!, df′! = JopCircShift_df′!, dom = spc, rng = spc, s = (shifts_forward=shifts,shifts_adjoint=map(-, shifts)))
end
export JopCircShift

JopCircShift_df!(d, m; shifts_forward, kwargs...) = circshift!(d, m, shifts_forward)
JopCircShift_df′!(m, d; shifts_adjoint, kwargs...) = circshift!(m, d, shifts_adjoint)
