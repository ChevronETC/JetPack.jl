
Jets.jl operator pack


<a id='Operator-documentation-1'></a>

# Operator documentation

- [`JetPack.JopAtan`](README.md#JetPack.JopAtan)
- [`JetPack.JopBlend`](README.md#JetPack.JopBlend)
- [`JetPack.JopDiagonal`](README.md#JetPack.JopDiagonal)
- [`JetPack.JopExp`](README.md#JetPack.JopExp)
- [`JetPack.JopGradient`](README.md#JetPack.JopGradient)
- [`JetPack.JopHighpass`](README.md#JetPack.JopHighpass)
- [`JetPack.JopImag`](README.md#JetPack.JopImag)
- [`JetPack.JopInterp`](README.md#JetPack.JopInterp)
- [`JetPack.JopLMO`](README.md#JetPack.JopLMO)
- [`JetPack.JopLaplacian`](README.md#JetPack.JopLaplacian)
- [`JetPack.JopLog`](README.md#JetPack.JopLog)
- [`JetPack.JopNim`](README.md#JetPack.JopNim)
- [`JetPack.JopNormalize`](README.md#JetPack.JopNormalize)
- [`JetPack.JopPad`](README.md#JetPack.JopPad)
- [`JetPack.JopPermute`](README.md#JetPack.JopPermute)
- [`JetPack.JopPermutedims`](README.md#JetPack.JopPermutedims)
- [`JetPack.JopPow`](README.md#JetPack.JopPow)
- [`JetPack.JopReal`](README.md#JetPack.JopReal)
- [`JetPack.JopReghost`](README.md#JetPack.JopReghost)
- [`JetPack.JopShift`](README.md#JetPack.JopShift)
- [`JetPack.JopTaper`](README.md#JetPack.JopTaper)

<a id='JetPack.JopAtan' href='#JetPack.JopAtan'>#</a>
**`JetPack.JopAtan`** &mdash; *Function*.



```julia
F = JopAtan(spc, c)
```

where `F` is the shifted arc tangent operator F = arctan( x(t)/c ) + π/2 

x(t) = 0 corresponds to F = π/2 if |x(t)| ≤ c, then π/4 ≤ F ≤ 3π/4

<a id='JetPack.JopBlend' href='#JetPack.JopBlend'>#</a>
**`JetPack.JopBlend`** &mdash; *Function*.



```julia
A = JotOpBlend(T, nsamples, shottimes[, nsamples_blended])
```

where `A` is a shot mixing operation, so that if `d=Am`, then `m::Array{T,2}` is a common receiver gather with each trace corresponding to a different shot.  The excitation time of each are the entries of `shottimes::Array{Int,1}` in units of time samples.  The size of `m` is `(nsamples,length(shottimes))`.  The size of `d` is either `nsamples_blended` (if provided) or `maximum(shot_times)+nsamples-1`.

**Examples**

**blend two receiver gather traces:**

```julia
A = JotOpBlend(Float64, 128, [1,64])
m = rand(domain(A))
d = A*m # receiver gather with blended shots
```


!!! warning "Missing docstring."
    Missing docstring for `JopCircShift`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `JopDerivative`. Check Documenter's build log for details.


<a id='JetPack.JopDiagonal' href='#JetPack.JopDiagonal'>#</a>
**`JetPack.JopDiagonal`** &mdash; *Function*.



A = JopDiagonal(spc, d) A = JopDiagonal(d)

where `spc::JetSpace` is the domain/range of `A`, and `d::Array` or `d::Number` is the diagonal.  If `d<:Number`, then the diagonal of the matrix is constant and one must specify `spc`.

**Examples**

```julia
A = JopDiagonal(spc, d)
```

where `spc::JetSpace` is the domain and range of `A`, and `d<:Array` or `d<:Number`

```julia
A = JotOpDiagonal(d)
```

where `d<:AbstractArray`.  The domain and range of `A` are determined by the size and type of `d`.

```julia
A = JotOpDiagonal([1.0, 2.0, 3.0])
m = ones(domain(A))
d = A*m # d = [1.0 ; 2.0 ; 3.0]
```


!!! warning "Missing docstring."
    Missing docstring for `JopDifference`. Check Documenter's build log for details.


<a id='JetPack.JopExp' href='#JetPack.JopExp'>#</a>
**`JetPack.JopExp`** &mdash; *Function*.



```julia
F = JopExp(spc, c)
```

where `F` is the exponential operator e^(c/x) with domain and range given by `spc::JetSpace`, and scalar value c.

<a id='JetPack.JopGradient' href='#JetPack.JopGradient'>#</a>
**`JetPack.JopGradient`** &mdash; *Function*.



```julia
A = JopGradient(dom, δ)
```

Gradient for a 2D or 3D arrays. The range will have one more dimension than the domain and contain the components of the gradient in each dimension.

  * `dom::JetSpace{T,N}` is the domain of the operator.
  * `δ::NTuple{T,N}` is the grid spacing in each dimension.

**Examples:**

**2D**

```julia
nz,nx = 11,12
dom = JetSpace(Float32, nz, nx)
A = JopGradient(dom, (1.0,1.0))
size(domain(A)) # (11,12)
size(range(A))  # (11,12,2)
```

**3D**

```julia
nz,ny,nx = 11,12,13
dom = JetSpace(Float32, nz, ny, nx)
A = JopGradient(dom,(1.0,1.0,1.0))
size(domain(A)) # (11,12,13)
size(range(A))  # (11,12,13,3)
```

<a id='JetPack.JopHighpass' href='#JetPack.JopHighpass'>#</a>
**`JetPack.JopHighpass`** &mdash; *Function*.



```julia
A = JopHighpass(sp)
```

Build a 2D Highpass operator, which first apply a center weighted three point smoothing operator [0.25, 0.5, 0.25] to each dimension with a given iterations number (nit[1:2]) to get the low frequency component of the data, then subtract it from the original data to return the high pass component.  It works with `sp::JetSpace`, `nit::Array{Integer, 1}` and `A::JopHighpass`.  It can specify different number of iterations (nit) in x and z direction. For a seismic image, more iterations in z direction is often applied to obtain more vertical resolution.

For example,

```julia
A = JopHighpass(JetSpace(Float64,128,256), nit=(3,1)])
m = rand(domain(A))
d = A*m
```

<a id='JetPack.JopImag' href='#JetPack.JopImag'>#</a>
**`JetPack.JopImag`** &mdash; *Function*.



```julia
op = JopImag(dom)
```

Extract the imaginary part of a complex input array where `dom::JetSpace{<:Complex}`.

<a id='JetPack.JopInterp' href='#JetPack.JopInterp'>#</a>
**`JetPack.JopInterp`** &mdash; *Function*.



```julia
JopInterp(dom, rng)
```

Performs linear interpolation from dom::JetSpace to rng::JetSpace. It is often used to reduce dimensionality in FWI.  JopInterp is ported from the CVX frequency domain FWI tools in SeisSpace.

**Notes**

  * It is required that domain and range have the same dimensionality.
  * It is required that domain and range both have more than 2 points per dimension
  * It is required that the domain is coarser than the range.
  * It is assumed that domain and range have the same boundary (e.g. in 2D: xmin,xmax,zmin,zmax).

**Examples:**

**2D**

```julia
dom = JetSpace(Float32, 11, 11)
rng = JetSpace(Float32, 23, 23)
A = JopInterp(dom, rng)
y = A * rand(dom)
x = A' * rand(rng)
```

<a id='JetPack.JopLaplacian' href='#JetPack.JopLaplacian'>#</a>
**`JetPack.JopLaplacian`** &mdash; *Function*.



```julia
A = JopLaplacian(sp)
```

Build a 2D (d^2/dx^2+d^2/dy^2) or 3D (d^2/dx^2+d^2/dy^2+d^2/dz^2) Laplacian operator with three point centered finite difference stencil, with `sp::JetSpace` and `A::JopLaplacian`.

For example,

```julia
A = JopLaplacian(JetSpace(Float64,128,256))
m = rand(domain(A))
d = A*m
```


!!! warning "Missing docstring."
    Missing docstring for `JopLeakyIntegration`. Check Documenter's build log for details.


<a id='JetPack.JopLMO' href='#JetPack.JopLMO'>#</a>
**`JetPack.JopLMO`** &mdash; *Function*.



```julia
A = JopLMO(spc [;v=2.0, dx=1.0, dt=1.0])
```

2D linear moveout operator for domain and range `sp::JetSpace`.  The forward operator maps from flat events to dipping events using the parameters `v,dx,dt`, and:

```
(z_moveout) = (z + (dx/dt)/v * x)
```

**Notes**

  * It should be trivial to add a 3D implementation for this operator.

<a id='JetPack.JopLog' href='#JetPack.JopLog'>#</a>
**`JetPack.JopLog`** &mdash; *Function*.



```julia
F = JopLog(spc)
```

where `F` is the log operator with domain and range given by `spc::JetSpace`.

<a id='JetPack.JopNim' href='#JetPack.JopNim'>#</a>
**`JetPack.JopNim`** &mdash; *Function*.



```julia
F = JopNim(spc, c)
```

where `F` is the 'normalized integral method' operator F = int*0^t [ x(t) dt ] / int*0_T [ x(t) dt ]

<a id='JetPack.JopNormalize' href='#JetPack.JopNormalize'>#</a>
**`JetPack.JopNormalize`** &mdash; *Function*.



```julia
F = JopNormalize(spc, c)
```

where `F` is the 'normalize by maximuam value of integral' operator F(t) = x(t) / int_0^T [ x(t) dt ] Note the integration is along the 1st (fastest) dimension.

<a id='JetPack.JopPad' href='#JetPack.JopPad'>#</a>
**`JetPack.JopPad`** &mdash; *Function*.



```julia
A = JopPad(dom, pad..., [extend=false, accumulate=false])
```

where `dom::JetSpace` is the domain of `A`, and `pad::UnitRange...` determines the range of `A`. If `extend=false`, then the padded region is set to zero.  If `extend=true`, then the padded region is set from the boundary of the domain.  The `accumulate=true` option is specific to the `Ginsu` operation in JetPackWave, and should not be used unless you really know what you are doing.

**Examples:**

**1D**

```julia
A = JopPad(JetSpace(Float64,2), -1:3)
m = [1.0, 2.0]
d = A*m # d = [0.0, 0.0, 1.0, 2.0, 0.0]
A = JopPad(JetSpace(Float64,2), -1:3, extend=true)
d = A*m # d = [1.0, 1.0, 1.0, 2.0, 2.0]
```

**2D**

```julia
A = JopPad(JetSpace(Float64,2,2), -1:3, 1:3)
m = [11. 12. ; 21. 22.]
d = A*m # d = [0. 0. 11. 12. 0. ; 0. 0. 21. 22. 0. ; 0. 0. 0. 0. 0.]
```

**Notes:**

  * This operator may also be used for truncation
  * There may be overlap between the functionality of `JopPad` and `JopRestriction`, and it may be worth

thinking of how to consolidate them.

<a id='JetPack.JopPermute' href='#JetPack.JopPermute'>#</a>
**`JetPack.JopPermute`** &mdash; *Function*.



```julia
A = JopPermute(sp, dims, perm)
```

where `sp::JetSpace` is the domain and range of the operator, `dims::NTuple` is a tuple of dimensions that one wishes to permute, and `perm` is a tuple of arrays with permutation indices.

**Example 1:**

```julia
A = JopPermute(JetSpace(Float64,4,2),(1,),([3;1;2;4],))
m = [11 12 ; 21 22 ; 31 32 ; 41 42]
d = A*m # d=[31 32 ; 11 12 ; 21 22 ; 41 42]
```

**Example 2:**

```julia
A = JopPermute(JetSpace(Float64,3,2),(1,2),([3;2;1],[2;1]))
m = [11 12 ; 21 22 ; 31 32]
d = A*m # d = [32 31 ; 22 21 ; 12 11]
```

**Notes:**

Currently, this is only implmented for 1D, 2D and 3D arrays.

<a id='JetPack.JopPermutedims' href='#JetPack.JopPermutedims'>#</a>
**`JetPack.JopPermutedims`** &mdash; *Function*.



```julia
A = JopPermutedims(sp, perm)
```

where `sp::JetSpace` is the domain of the operator, and `perm` is an array specifying the permutation of array dimensions.  The range of the operator is inferred from `sp` and `perm`.

**Example:**

```julia
A = JopPermutedims(JetSpace(Float64,3,2),[2;1])
m = [1 2 3 ; 4 5 6]
d = A*m # d = [1 4 ; 2 5 ; 3 6]
```

<a id='JetPack.JopPow' href='#JetPack.JopPow'>#</a>
**`JetPack.JopPow`** &mdash; *Function*.



```julia
F = JopPow(spc, c, a)
```

where `F` is the power operator (x/c)^a with domain and range given by `spc::JetSpace`, and scalar values c, a.


!!! warning "Missing docstring."
    Missing docstring for `JopProjection`. Check Documenter's build log for details.


<a id='JetPack.JopReal' href='#JetPack.JopReal'>#</a>
**`JetPack.JopReal`** &mdash; *Function*.



```julia
op = JopReal(dom)
```

Extract the real part of a complex input array where `dom::JetSpace{<:Complex}`.

<a id='JetPack.JopReghost' href='#JetPack.JopReghost'>#</a>
**`JetPack.JopReghost`** &mdash; *Function*.



```julia
A = JopReghost(spc, zt, zm, dt, dx)
```

Build a 2D receiver-side reghosting and redatuming operator in the FK domain using simple phase-shift operations (e.g. Posthumus, 1993).  `spc::JetSpaceSymmetric` is the domain and range of `A`. `zt` is target depth where the receivers will be re-dataumed to, and `zm` is the measurement depth where the recording is made. `dt` is the time sampling interval and `dx` is the receiver sampling interval.

Since this operator is built in the **FK** domain, `spc` has symmetry in the frequency dimension.  Hence `spc::JotSpaceSymmetric`.  We provide a convenience method:

```
sp = symspace(JopReghost, T, nw, nk)
```

where `T` is either `Float32` or `Float64`, `nw` are the number of frequency samples, and `nk` are the number of wavenumber samples.

**Example:**

```julia
P = JopPad(JetSpace(Float64,nt,nx), 1:2nt, 1:nx)
F = JopFft(range(P))
G = JopReghost(range(F), 20.0, 20.0, .004, 10.0)
A = P'*F'*G*F*P
data_reghost = A * data_noghost
```

Note that in the above example,

```julia
range(F) == symspace(JopReghost_df!, Float64, size(range(P))...)
```


!!! warning "Missing docstring."
    Missing docstring for `JopReoughness`. Check Documenter's build log for details.


<a id='JetPack.JopShift' href='#JetPack.JopShift'>#</a>
**`JetPack.JopShift`** &mdash; *Function*.



A = JopShift(spc, b) A = JopShift(b)

A(x) = x + b

where `spc::JetSpace` is the domain/range of `A`, and `b::Array` or `b::Number` is the affine translation or shift.  If `b<:Number`, then shift is constant and one must specify `spc`.

**Examples**

```julia
A = JopShift(spc, b)
```

where `spc::JetSpace` is the domain and range of `A`, and `b<:Array` or `b<:Number`

```julia
A = JopShift(b)
```

where `b<:AbstractArray`.  The domain and range of `A` are determined by the size and type of `b`.

```julia
A = JopShift([1.0, 2.0, 3.0])
m = ones(domain(A))
d = A*m # d = [2.0 ; 3.0 ; 4.0]
```

<a id='JetPack.JopTaper' href='#JetPack.JopTaper'>#</a>
**`JetPack.JopTaper`** &mdash; *Function*.



```julia
A = JopTaper(spc, dims, frac[, frac_end; mode=:normal, taper=(:cosine,:cosine)])
```

The linear operator `A` tapers the edges of all or some subset of the dimensions of an array that belongs to `spc::JetSpace`.  The dimensions that are tapered are given by `dims`, an `Int` tuple. `A` will taper the beginning and/or end edges of each specified array dimension. The size of the beginning and end tapers are determined as a fraction of the length of that array dimension, and this fraction is set using the tuple `frac`, and (optionally) `frac_end`.

If `frac_end::NTuple{M,Float64}` is not set, then the ith entry of `frac::NTuple{M,Float64}` determines the length of the center portion of the ith dimension that is not tapered.  If `frac_end::NTuple{M,NTuple{2,Float64}}` is set, then the length of the begining taper of the ith array dimension is determiend by `frac[i]`, and the length of the end taper is determined by `frac_end[i]`.

The optional named argument `mode` can be used if the taper is applied to a dimension that corresponds to FFT ordering where the edges are assumed to be at the center and left ends of the dimension.

The optional named argument `taper` is a tuple specifying what type of taper to use at each end.  Available tapers are :cosine and :heviside.

**Examples:**

**taper for a space containing 1D arrays**

```julia
A = JopTaper(JetSpace(Float64,10), (1,), (.75,))
m = ones(domain(A))
d = A*m
```

**taper for a space containing 1D arrays, and only taper at the end**

```julia
A = JopTaper(JetSpace(Float64,10), (1,), (0.0,), (0.25,))
m = ones(domain(A))
d = A*m
```

**taper for a space containing 2D arrays, where both dimensions are tapered**

```julia
A = JopTaper(JetSpace(Float64,10,11), (1,2), (0.75,0.5))
m = ones(domain(A))
d = A*m
```

**taper for a space containing 3D arrays, where two of the three dimensions are tapered**

```julia
A = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.75,0.5))
m = ones(domain(A))
d = A*m
```

In the next example, array dimension `1` is tapered at the end, and array dimension `2` is tapered at the beginning and end.

```julia
A = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.0,0.25), (0.25,0.25))
m = ones(domain(A))
d = A*m
```


!!! warning "Missing docstring."
    Missing docstring for `JopTranslation`. Check Documenter's build log for details.


