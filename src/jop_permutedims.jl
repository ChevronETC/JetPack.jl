"""
    A = JopPermutedims(sp, perm)

where `sp::JetSpace` is the domain of the operator, and `perm` is an array
specifying the permutation of array dimensions.  The range of the operator
is inferred from `sp` and `perm`.

# Example:
```julia
A = JopPermutedims(JetSpace(Float64,3,2),[2;1])
m = [1 2 3 ; 4 5 6]
d = A*m # d = [1 4 ; 2 5 ; 3 6]
```
"""
function JopPermutedims(dom, perm)
    shpdom = size(dom)
    shprng = ntuple(i->shpdom[perm[i]], length(shpdom))
    rng = JetSpace(eltype(dom), shprng)

    P = zeros(Int,length(perm),length(perm))
    for i = 1:length(perm)
        P[i,perm[i]] = 1
    end
    iperm = P'*[1:length(shpdom);]

    JopLn(dom = dom, rng = rng, df! =  JopPermutedims_df!, df′! = JopPermutedims_df′!, s = (perm=perm, iperm=iperm))
end

export JopPermutedims

JopPermutedims_df!(d, m; perm, kwargs...) = permutedims!(d, m, perm)
JopPermutedims_df′!(m, d; iperm, kwargs...) = permutedims!(m, d, iperm)
