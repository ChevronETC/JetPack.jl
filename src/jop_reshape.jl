"""
    A = JopReshape(dom, rng)

Reshape an array that belongs in `dom::JetAbstractSpace` to one that belongs to `rng::JetAbstractSpace`

# Example
```julia
using Jets, JetPack
dom = JetSpace(Float32,10,20)
rng = JetSpace(Float32,200)
A = JopReshape(JetSpace(Float32,10,20), JetSpace(Float32,200))
x = rand(domain(A))
y = A*x
```
"""
JopReshape(dom::JetAbstractSpace, rng::JetAbstractSpace) = JopLn(dom = dom, rng = rng, df! = JopReshape_df!)

export JopReshape

function JopReshape_df!(d, m; kwargs...)
    d[:] .= m[:]
    d
end
