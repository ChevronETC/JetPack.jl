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
function JopAffine(shift::AbstractArray{T}) where {T}
    spc = JetSpace(T, size(shift))
    JopNl(f! = JopAffine_f!, df! = JopAffine_df!, df′! = JopAffine_df′!, dom = spc, rng = spc, s = (shift=shift,))
end

JopAffine(spc::JetAbstractSpace, shift) = JopNl(f! = JopAffine_f!, df! = JopAffine_df!, df′! = JopAffine_df′!, dom = spc, rng = spc, s = (shift=shift,))
export JopAffine

JopAffine_f!(d, m; shift, kwargs...) = d .=  m .+ shift
JopAffine_df!(δd, δm; mₒ, shift, kwargs...) = δd .=  δm #Derivative is identity
JopAffine_df′!(δm, δd; mₒ, shift, kwargs...) = δm .=  δd #Adjoint of derivative is identity
