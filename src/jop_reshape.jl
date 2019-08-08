JopReshape(dom::JetAbstractSpace, rng::JetAbstractSpace) = JopLn(dom = dom, rng = rng, df! = JopReshape_df!)

export JopReshape

function JopReshape_df!(d, m; kwargs...)
    d[:] .= m[:]
    d
end
