JopReshape(dom::JetAbstractSpace, rng::JetAbstractSpace) = JopLn(dom = dom, rng = rng, df! = JopReshape_df!, df′! = JopReshape_df′!)

export JopReshape

JopReshape_df!(d, m; kwargs...) = d[:] .= m[:]
JopReshape_df′!(m, d; kwargs...) = m[:] .= d[:]
