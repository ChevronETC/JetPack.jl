"""
A = JopShift(spc, b)
A = JopShift(b)

A(x) = x + b

where `spc::JetSpace` is the domain/range of `A`, and `b::Array` or `b::Number`
is the affine translation or shift.  If `b<:Number`, then shift is constant
and one must specify `spc`.

# Examples
```julia
A = JopShift(spc, b)
```
where `spc::JetSpace` is the domain and range of `A`, and `b<:Array` or
`b<:Number`

```julia
A = JopShift(b)
```
where `b<:AbstractArray`.  The domain and range of `A` are determined by the
size and type of `b`.

```julia
A = JopShift([1.0, 2.0, 3.0])
m = ones(domain(A))
d = A*m # d = [2.0 ; 3.0 ; 4.0]
```
"""
JopShift(spc::JetAbstractSpace, shift) = JopNl(f! = JopShift_f!, df! = JopShift_df!, dom = spc, rng = spc, s = (shift=shift,))
JopShift(shift::AbstractArray{T}) where {T} = JopShift( JetSpace(T, size(shift)),shift)
export JopShift

JopShift_f!(d, m; shift, kwargs...) = d .=  m .+ shift
JopShift_df!(δd, δm; mₒ, shift, kwargs...) = δd .=  δm #Derivative is identity
