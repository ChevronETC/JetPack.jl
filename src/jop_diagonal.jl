"""
A = JopDiagonal(spc, d)
A = JopDiagonal(d)

where `spc::JetSpace` is the domain/range of `A`, and `d::Array` or `d::Number`
is the diagonal.  If `d<:Number`, then the diagonal of the matrix is constant
and one must specify `spc`.

# Examples
```julia
A = JopDiagonal(spc, d)
```
where `spc::JetSpace` is the domain and range of `A`, and `d<:Array` or
`d<:Number`

```julia
A = JotOpDiagonal(d)
```
where `d<:AbstractArray`.  The domain and range of `A` are determined by the
size and type of `d`.

```julia
A = JotOpDiagonal([1.0, 2.0, 3.0])
m = ones(domain(A))
d = A*m # d = [1.0 ; 2.0 ; 3.0]
```
"""
function JopDiagonal(diagonal::AbstractArray)
    spc = space(diagonal)
    JopLn(df! = JopDiagonal_df!, df′! = JopDiagonal_df′!, dom = spc, rng = spc, s = (diagonal=diagonal,))
end
JopDiagonal(spc::JetAbstractSpace, diagonal) = JopLn(df! = JopDiagonal_df!, df′! = JopDiagonal_df′!, dom = spc, rng = spc, s = (diagonal=diagonal,))
export JopDiagonal

JopDiagonal_df!(d, m; diagonal, kwargs...) = d .= diagonal .* m
JopDiagonal_df′!(m, d; diagonal, kwargs...) = m .= conj(diagonal) .* d
