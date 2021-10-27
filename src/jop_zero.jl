"""
    A = JopZeros(dom::JetAbstractSpace,rng::JetAbstractSpace)

`A*m` transforms m ∈ dom to the cannonical zero vector ∈ rng    
"""
function JopZero(dom::JetAbstractSpace, rng::JetAbstractSpace) 
    JopLn(dom = dom, rng = rng, df! = JopZero_df!, df′! = JopZero_df′!, s=(dom=dom, rng=rng))
end

JopZero_df!(d, m; rng::JetAbstractSpace, kwargs...) = zeros(rng)
JopZero_df′!(m, d; dom::JetAbstractSpace, kwargs...) = zeros(dom)

export JopZero

