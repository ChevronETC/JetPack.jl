# JetPack

A somewhat hap-hazard set of operators for [Jets.jl](https://github.com/ChevronETC/Jets.jl).

## Example
```julia
using Jets, JetPack, Random
R = JetSpace(Float64,128,128)
x = ones(R)
A = JopRestriction(R, randperm(R, 512))
y = A*x
z = A'*y
```
One can now plot `x` and `z` to see the effect of the restriction.  For example,
using PyPlot:
```julia
using PyPlot
figure(1);clf()
subplot(121);imshow(x)
subplot(122);imshow(z) 
```