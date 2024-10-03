"""
    A = JopSpreadStack(dom, rng)

where `dom::JetSpace{T,1}` is the 1D domain of `A`, and `rng::JetSpace{T,2}` is the 2D range of 
`A`, with the range having one higher dimension than the domain. The forward operator copies 
(spreads) the domain vector across the range vector, and the adjoint operator increments (stacks)
across the range vector into the domain vector.

# Examples:

## 1D
```julia
A = JopSpreadStack(JetSpace(T,11), JetSpace(T,11,5))
x = rand(domain(A))
y = A * x
z = A' * y
```
"""
function JopSpreadStack(dom::JetSpace{T,1}, rng::JetSpace{T,2}) where {T}
    JopLn(dom = dom, rng = rng, df! = JopSpreadStack_df!, df′! = JopSpreadStack_df′!)
end
export JopSpreadStack

function JopSpreadStack_df!(d::AbstractArray{T,2}, m::AbstractArray{T,1}; kwargs...) where {T}
    d .= 0
    for k ∈ 1:size(d,2)
        d[:,k] .= m[:]
    end
    d
end

function JopSpreadStack_df′!(m::AbstractArray{T,1}, d::AbstractArray{T,2}; kwargs...) where {T}
    m .= 0
    for k ∈ 1:size(d,2)
        m[:] .+= d[:,k]
    end
    m
end
