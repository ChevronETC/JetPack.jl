function JopCircShift(spc, shifts)
    df!(d, m; shifts_forward, kwargs...) = circshift!(d, m, shifts_forward)
    df′!(m, d; shifts_adjoint, kwargs...) = circshift!(m, d, shifts_adjoint)

    JopLn(df! = df!, df′! = df′!, dom = spc, rng = spc, s = (shifts_forward=shifts,shifts_adjoint=map(-, shifts)))
end
export JopCircShift
