"""
A = JopAffine(spc, b)
A = JopAffine(b)

A(x) = x + b

where `spc::JetSpace` is the domain/range of `A`, and `b::Array` or `b::Number`
is the affine translation or shift.  If `b<:Number`, then shift is constant
and one must specify `spc`.

# Examples
```julia
A = JopAffine(spc, b)
```
where `spc::JetSpace` is the domain and range of `A`, and `b<:Array` or
`b<:Number`

```julia
A = JotOpDiagonal(b)
```
where `b<:AbstractArray`.  The domain and range of `A` are determined by the
size and type of `b`.

```julia
A = JopAffine([1.0, 2.0, 3.0])
m = ones(domain(A))
d = A*m # d = [2.0 ; 3.0 ; 4.0]
```
"""
JopAffine(spc::JetAbstractSpace, shift) = JopNl(f! = JopAffine_f!, df! = JopAffine_df!, dom = spc, rng = spc, s = (shift=shift,))
JopAffine(shift::AbstractArray{T}) where {T} = JopAffine( JetSpace(T, size(shift)),shift)
export JopAffine

JopAffine_f!(d, m; shift, kwargs...) = d .=  m .+ shift
JopAffine_df!(δd, δm; mₒ, shift, kwargs...) = δd .=  δm #Derivative is identity
