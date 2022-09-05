var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [JetPack]\nOrder = [:function]","category":"page"},{"location":"reference/#JetPack.JopAtan-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}, Tuple{Jets.JetSpace{T}, Any}} where T","page":"Reference","title":"JetPack.JopAtan","text":"F = JopAtan(spc, c)\n\nwhere F is the shifted arc tangent operator F = arctan( x(t)/c ) + π/2 \n\nx(t) = 0 corresponds to F = π/2 if |x(t)| ≤ c, then π/4 ≤ F ≤ 3π/4\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopBlend","page":"Reference","title":"JetPack.JopBlend","text":"A = JotOpBlend(T, nsamples, shottimes[, nsamples_blended])\n\nwhere A is a shot mixing operation, so that if d=Am, then m::Array{T,2} is a common receiver gather with each trace corresponding to a different shot.  The excitation time of each are the entries of shottimes::Array{Int,1} in units of time samples.  The size of m is (nsamples,length(shottimes)).  The size of d is either nsamples_blended (if provided) or maximum(shot_times)+nsamples-1.\n\nExamples\n\nblend two receiver gather traces:\n\nA = JotOpBlend(Float64, 128, [1,64])\nm = rand(domain(A))\nd = A*m # receiver gather with blended shots\n\n\n\n\n\n","category":"function"},{"location":"reference/#JetPack.JopCircShift-Tuple{Any, Any}","page":"Reference","title":"JetPack.JopCircShift","text":"A = JopCircShift(R::JetSpace, shifts)\n\nd=A*m is equivalent to circshift(d, m, shifts), and where m and d are in the domain/range of A.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopDerivative-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}} where T","page":"Reference","title":"JetPack.JopDerivative","text":"A = JopDerivative(R::JetSpace[; dim=1, accuracy=4, delta=1.0)\n\nA*m is the centered finite different approximation to the derivative of m along dimension dim.  The accuracy of the approximation is either 4 for a 4th order accurate estimate of the derivative, or 8 for an 8th order accurate derivative.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopDiagonal-Tuple{AbstractArray}","page":"Reference","title":"JetPack.JopDiagonal","text":"A = JopDiagonal(spc, d) A = JopDiagonal(d)\n\nwhere spc::JetSpace is the domain/range of A, and d::Array or d::Number is the diagonal.  If d<:Number, then the diagonal of the matrix is constant and one must specify spc.\n\nExamples\n\nA = JopDiagonal(spc, d)\n\nwhere spc::JetSpace is the domain and range of A, and d<:Array or d<:Number\n\nA = JotOpDiagonal(d)\n\nwhere d<:AbstractArray.  The domain and range of A are determined by the size and type of d.\n\nA = JotOpDiagonal([1.0, 2.0, 3.0])\nm = ones(domain(A))\nd = A*m # d = [1.0 ; 2.0 ; 3.0]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopDifference","page":"Reference","title":"JetPack.JopDifference","text":"A = JopDifference(R::JetSpace[, dim=1])\n\nA*m is similar to diff(a,dims=dim) where a is a vector in the space R. In other words, it is the one-sided difference of a along the dimension dim.\n\n\n\n\n\n","category":"function"},{"location":"reference/#JetPack.JopErf-Tuple{Jets.JetSpace}","page":"Reference","title":"JetPack.JopErf","text":"F = JopErf(spc)\n\nwhere F is the error function operator with domain and range given by spc::JetSpace. we use 'erf'(z) = 2/√π ∫_{0}^{z} \u001bxp{-t^2}dt, the derivative of which is 2/√π  \u001bxp{-z^2}. We expect the domain and range to be real.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopExp-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}, Tuple{Jets.JetSpace{T}, Any}} where T","page":"Reference","title":"JetPack.JopExp","text":"F = JopExp(spc, c)\n\nwhere F is the exponential operator e^(c/x) with domain and range given by spc::JetSpace, and scalar value c.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopGradient-Union{Tuple{N}, Tuple{T}, Tuple{Jets.JetSpace{T, N}, Tuple{Vararg{T, N}}}} where {T, N}","page":"Reference","title":"JetPack.JopGradient","text":"A = JopGradient(dom, δ)\n\nGradient for a 2D or 3D arrays. The range will have one more dimension than the domain and contain the components of the gradient in each dimension.\n\ndom::JetSpace{T,N} is the domain of the operator.\nδ::NTuple{T,N} is the grid spacing in each dimension.\n\nExamples:\n\n2D\n\nnz,nx = 11,12\ndom = JetSpace(Float32, nz, nx)\nA = JopGradient(dom, (1.0,1.0))\nsize(domain(A)) # (11,12)\nsize(range(A))  # (11,12,2)\n\n3D\n\nnz,ny,nx = 11,12,13\ndom = JetSpace(Float32, nz, ny, nx)\nA = JopGradient(dom,(1.0,1.0,1.0))\nsize(domain(A)) # (11,12,13)\nsize(range(A))  # (11,12,13,3)\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopHighpass-Tuple{Jets.JetSpace}","page":"Reference","title":"JetPack.JopHighpass","text":"A = JopHighpass(sp)\n\nBuild a 2D Highpass operator, which first apply a center weighted three point smoothing operator [0.25, 0.5, 0.25] to each dimension with a given iterations number (nit[1:2]) to get the low frequency component of the data, then subtract it from the original data to return the high pass component.  It works with sp::JetSpace, nit::Array{Integer, 1} and A::JopHighpass.  It can specify different number of iterations (nit) in x and z direction. For a seismic image, more iterations in z direction is often applied to obtain more vertical resolution.\n\nFor example,\n\nA = JopHighpass(JetSpace(Float64,128,256), nit=(3,1)])\nm = rand(domain(A))\nd = A*m\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopImag-Union{Tuple{Jets.JetSpace{Complex{T}}}, Tuple{T}} where T","page":"Reference","title":"JetPack.JopImag","text":"op = JopImag(dom)\n\nExtract the imaginary part of a complex input array where dom::JetSpace{<:Complex}.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopInterp-Tuple{Jets.JetAbstractSpace, Jets.JetAbstractSpace}","page":"Reference","title":"JetPack.JopInterp","text":"JopInterp(dom, rng)\n\nPerforms linear interpolation from dom::JetSpace to rng::JetSpace. It is often used to reduce dimensionality in FWI.  JopInterp is ported from the CVX frequency domain FWI tools in SeisSpace.\n\nNotes\n\nIt is required that domain and range have the same dimensionality.\nIt is required that domain and range both have more than 2 points per dimension\nIt is required that the domain is coarser than the range.\nIt is assumed that domain and range have the same boundary (e.g. in 2D: xmin,xmax,zmin,zmax).\n\nExamples:\n\n2D\n\ndom = JetSpace(Float32, 11, 11)\nrng = JetSpace(Float32, 23, 23)\nA = JopInterp(dom, rng)\ny = A * rand(dom)\nx = A' * rand(rng)\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopLMO-Union{Tuple{T}, Tuple{Jets.JetSpace{T}, Any, Any, Any}} where T","page":"Reference","title":"JetPack.JopLMO","text":"A = JopLMO(spc [;v=2.0, dx=1.0, dt=1.0])\n\n2D linear moveout operator for domain and range sp::JetSpace.  The forward operator maps from flat events to dipping events using the parameters v,dx,dt, and:\n\n(z_moveout) = (z + (dx/dt)/v * x)\n\nNotes\n\nIt should be trivial to add a 3D implementation for this operator.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopLaplacian-Tuple{Jets.JetSpace}","page":"Reference","title":"JetPack.JopLaplacian","text":"A = JopLaplacian(sp)\n\nBuild a 2D (d^2/dx^2+d^2/dy^2) or 3D (d^2/dx^2+d^2/dy^2+d^2/dz^2) Laplacian operator with three point centered finite difference stencil, with sp::JetSpace and A::JopLaplacian.\n\nFor example,\n\nA = JopLaplacian(JetSpace(Float64,128,256))\nm = rand(domain(A))\nd = A*m\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopLog-Tuple{Jets.JetSpace}","page":"Reference","title":"JetPack.JopLog","text":"F = JopLog(spc)\n\nwhere F is the log operator with domain and range given by spc::JetSpace.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopMix-Union{Tuple{N}, Tuple{T}, Tuple{Jets.JetAbstractSpace{T, N}, Tuple{Vararg{Int64, N}}}} where {T, N}","page":"Reference","title":"JetPack.JopMix","text":"JopMix(s::JetAbstractSpace{T,N}, nmix::NTuple{N,Int})\n\n2D spatial mix, and is used to add smoothness to a domain.  JopMix is ported from the CVX frequency domain FWI tools in SeisSpace.\n\nNotes:\n\nThis appears to over-lap in funtionality with JopRoughness.\n\nExamples:\n\n2D\n\nnz, nx = 11,11\nspc = JetSpace(Float32, nz, nx)\nA = JopMix(spc, (5,5))\ny = A * rand(dom)\nx = A' * rand(rng)\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopNim-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}} where T","page":"Reference","title":"JetPack.JopNim","text":"F = JopNim(spc, c)\n\nwhere F is the 'normalized integral method' operator F = int0^t [ x(t) dt ] / int0_T [ x(t) dt ]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopNormalize-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}} where T","page":"Reference","title":"JetPack.JopNormalize","text":"F = JopNormalize(spc, c)\n\nwhere F is the 'normalize by maximuam value of integral' operator F(t) = x(t) / int_0^T [ x(t) dt ] Note the integration is along the 1st (fastest) dimension.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopPad-Union{Tuple{T}, Tuple{Jets.JetSpace{T}, Vararg{UnitRange}}} where T","page":"Reference","title":"JetPack.JopPad","text":"A = JopPad(dom, pad..., [extend=false, accumulate=false])\n\nwhere dom::JetSpace is the domain of A, and pad::UnitRange... determines the range of A. If extend=false, then the padded region is set to zero.  If extend=true, then the padded region is set from the boundary of the domain.  The accumulate=true option is specific to the Ginsu operation in JetPackWave, and should not be used unless you really know what you are doing.\n\nExamples:\n\n1D\n\nA = JopPad(JetSpace(Float64,2), -1:3)\nm = [1.0, 2.0]\nd = A*m # d = [0.0, 0.0, 1.0, 2.0, 0.0]\nA = JopPad(JetSpace(Float64,2), -1:3, extend=true)\nd = A*m # d = [1.0, 1.0, 1.0, 2.0, 2.0]\n\n2D\n\nA = JopPad(JetSpace(Float64,2,2), -1:3, 1:3)\nm = [11. 12. ; 21. 22.]\nd = A*m # d = [0. 0. 11. 12. 0. ; 0. 0. 21. 22. 0. ; 0. 0. 0. 0. 0.]\n\nNotes:\n\nThis operator may also be used for truncation\nThere may be overlap between the functionality of JopPad and JopRestriction, and it may be worth\n\nthinking of how to consolidate them.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopPermute-Tuple{Any, Any, Any}","page":"Reference","title":"JetPack.JopPermute","text":"A = JopPermute(sp, dims, perm)\n\nwhere sp::JetSpace is the domain and range of the operator, dims::NTuple is a tuple of dimensions that one wishes to permute, and perm is a tuple of arrays with permutation indices.\n\nExample 1:\n\nA = JopPermute(JetSpace(Float64,4,2),(1,),([3;1;2;4],))\nm = [11 12 ; 21 22 ; 31 32 ; 41 42]\nd = A*m # d=[31 32 ; 11 12 ; 21 22 ; 41 42]\n\nExample 2:\n\nA = JopPermute(JetSpace(Float64,3,2),(1,2),([3;2;1],[2;1]))\nm = [11 12 ; 21 22 ; 31 32]\nd = A*m # d = [32 31 ; 22 21 ; 12 11]\n\nNotes:\n\nCurrently, this is only implmented for 1D, 2D and 3D arrays.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopPermutedims-Tuple{Any, Any}","page":"Reference","title":"JetPack.JopPermutedims","text":"A = JopPermutedims(sp, perm)\n\nwhere sp::JetSpace is the domain of the operator, and perm is an array specifying the permutation of array dimensions.  The range of the operator is inferred from sp and perm.\n\nExample:\n\nA = JopPermutedims(JetSpace(Float64,3,2),[2;1])\nm = [1 2 3 ; 4 5 6]\nd = A*m # d = [1 4 ; 2 5 ; 3 6]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopPow-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}, Tuple{Jets.JetSpace{T}, Any}, Tuple{Jets.JetSpace{T}, Any, Any}} where T","page":"Reference","title":"JetPack.JopPow","text":"F = JopPow(spc, c, a)\n\nwhere F is the power operator (x/c)^a with domain and range given by spc::JetSpace, and scalar values c, a.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopProjection-Union{Tuple{AbstractArray{T, N}}, Tuple{N}, Tuple{T}} where {T, N}","page":"Reference","title":"JetPack.JopProjection","text":"A = JopProjection(u::AbstractArray)\n\nA*m is the projection of vector m onto the vector u.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopReal-Union{Tuple{Jets.JetSpace{Complex{T}}}, Tuple{T}} where T","page":"Reference","title":"JetPack.JopReal","text":"op = JopReal(dom)\n\nExtract the real part of a complex input array where dom::JetSpace{<:Complex}.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopReghost-Tuple{Jets.JetSSpace, Any, Any, Any, Any}","page":"Reference","title":"JetPack.JopReghost","text":"A = JopReghost(spc, zt, zm, dt, dx)\n\nBuild a 2D receiver-side reghosting and redatuming operator in the FK domain using simple phase-shift operations (e.g. Posthumus, 1993).  spc::JetSpaceSymmetric is the domain and range of A. zt is target depth where the receivers will be re-dataumed to, and zm is the measurement depth where the recording is made. dt is the time sampling interval and dx is the receiver sampling interval.\n\nSince this operator is built in the FK domain, spc has symmetry in the frequency dimension.  Hence spc::JotSpaceSymmetric.  We provide a convenience method:\n\nsp = symspace(JopReghost, T, nw, nk)\n\nwhere T is either Float32 or Float64, nw are the number of frequency samples, and nk are the number of wavenumber samples.\n\nExample:\n\nP = JopPad(JetSpace(Float64,nt,nx), 1:2nt, 1:nx)\nF = JopFft(range(P))\nG = JopReghost(range(F), 20.0, 20.0, .004, 10.0)\nA = P'*F'*G*F*P\ndata_reghost = A * data_noghost\n\nNote that in the above example,\n\nrange(F) == symspace(JopReghost_df!, Float64, size(range(P))...)\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopRemoveDC-Union{Tuple{Jets.JetAbstractSpace{T, N}}, Tuple{N}, Tuple{T}} where {T, N}","page":"Reference","title":"JetPack.JopRemoveDC","text":"JopRemoveDC(s::JetAbstractSpace{T,N})\n\nF(x) = x - zero frequency component of x along the fastest dimension Uses Fourier transform, as the DC component differs from the mean as a function of interval length and zero padding. ```\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopReshape-Tuple{Jets.JetAbstractSpace, Jets.JetAbstractSpace}","page":"Reference","title":"JetPack.JopReshape","text":"A = JopReshape(dom, rng)\n\nReshape an array that belongs in dom::JetAbstractSpace to one that belongs to rng::JetAbstractSpace\n\nExample\n\nusing Jets, JetPack\ndom = JetSpace(Float32,10,20)\nrng = JetSpace(Float32,200)\nA = JopReshape(JetSpace(Float32,10,20), JetSpace(Float32,200))\nx = rand(domain(A))\ny = A*x\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopRestriction-Union{Tuple{T}, Tuple{Jets.JetAbstractSpace{T}, Vector}} where T","page":"Reference","title":"JetPack.JopRestriction","text":"A = JopRestriction(dom[, rng], indices)\n\nApply a restriction operator mapping from dom::JetSpace to rng::JetSpace, and using a one-dimensional array of indices where indices is one of:\n\nindices::Vector{Int,1} - for restriction of 1-D arrays\nindices::Vector{NTuple{N,Int}} - for restriction of N-D arrays (N>1)\n\nNote that in the case where the domain is N-dimension, the range is always stored using a 1-dimension array.\n\nexample (1-D):\n\nA = JopRestriction(JotSpace(Float64,4),[1,3])\nm = [1.0;2.0;3.0;4.0]\nd = A*m # d=[1.0;3.0]\n\nexample (2-D):\n\nA = JopRestriction(JotSpace(Float64,2,2),[(1,2),(2,2)])\nm = [1.0 2.0;3.0 4.0]\nd = A*m # d=[2.0,4.0]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopRoughness-Tuple{Jets.JetSpace, Int64}","page":"Reference","title":"JetPack.JopRoughness","text":"A = JopRoughness(sp, dim, w)\n\nA  can be used to regularize an optimization problem, applying a penalty to models that are non-smooth.  This is similar to a finite difference operator, but here no care is taken to ensure that A computes a derivative.  dim is the dimension that one wishes to smooth and w is the half-width of the smoother. The width of the window determines how strong the penalty is for being non-smooth.\n\nFor example, the form of A for w=1 is,\n\nA=\n[\n(1/3-2/3) 1/3       0         0         ;\n1/3       (1/3-3/3) 1/3       0         ;\n0         1/3       (1/3-3/3) 1/3       ;\n0         0         1/3       (1/3-2/3)\n]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopShift-Tuple{Jets.JetAbstractSpace, Any}","page":"Reference","title":"JetPack.JopShift","text":"A = JopShift(spc, b) A = JopShift(b)\n\nA(x) = x + b\n\nwhere spc::JetSpace is the domain/range of A, and b::Array or b::Number is the affine translation or shift.  If b<:Number, then shift is constant and one must specify spc.\n\nExamples\n\nA = JopShift(spc, b)\n\nwhere spc::JetSpace is the domain and range of A, and b<:Array or b<:Number\n\nA = JopShift(b)\n\nwhere b<:AbstractArray.  The domain and range of A are determined by the size and type of b.\n\nA = JopShift([1.0, 2.0, 3.0])\nm = ones(domain(A))\nd = A*m # d = [2.0 ; 3.0 ; 4.0]\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopSigmoid-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}, Tuple{Jets.JetSpace{T}, Any}} where T","page":"Reference","title":"JetPack.JopSigmoid","text":"F = JopSigmoid(spc, c)\n\nwhere F is the sigmoid operator 1/(1+e^(-c/x)) with domain and range given by spc::JetSpace, and scalar value c.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopTanh-Union{Tuple{Jets.JetSpace{T}}, Tuple{T}, Tuple{Jets.JetSpace{T}, Any}} where T","page":"Reference","title":"JetPack.JopTanh","text":"F = JopTanh(spc)\n\nwhere F is the hyperbolic tangent operator with domain and range given by spc::JetSpace. we use 'tanh'(cx) = (\u001bxp{cx}-\u001bxp{-cx})/(\u001bxp{cx}+\u001bxp{-cx}), the derivative of which is 1-tanh(cx)^2. We expect the domain and range to be real.\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopTaper-Union{Tuple{M}, Tuple{Jets.JetAbstractSpace, Tuple{Vararg{Int64, M}}, Tuple{Vararg{Real, M}}, Tuple{Vararg{Real, M}}}} where M","page":"Reference","title":"JetPack.JopTaper","text":"A = JopTaper(spc, dims, frac[, frac_end; mode=:normal, taper=(:cosine,:cosine)])\n\nThe linear operator A tapers the edges of all or some subset of the dimensions of an array that belongs to spc::JetSpace.  The dimensions that are tapered are given by dims, an Int tuple. A will taper the beginning and/or end edges of each specified array dimension. The size of the beginning and end tapers are determined as a fraction of the length of that array dimension, and this fraction is set using the tuple frac, and (optionally) frac_end.\n\nIf frac_end::NTuple{M,Float64} is not set, then the ith entry of frac::NTuple{M,Float64} determines the length of the center portion of the ith dimension that is not tapered.  If frac_end::NTuple{M,NTuple{2,Float64}} is set, then the length of the begining taper of the ith array dimension is determiend by frac[i], and the length of the end taper is determined by frac_end[i].\n\nThe optional named argument mode can be used if the taper is applied to a dimension that corresponds to FFT ordering where the edges are assumed to be at the center and left ends of the dimension.\n\nThe optional named argument taper is a tuple specifying what type of taper to use at each end.  Available tapers are :cosine and :heaviside.\n\nExamples:\n\ntaper for a space containing 1D arrays\n\nA = JopTaper(JetSpace(Float64,10), (1,), (.75,))\nm = ones(domain(A))\nd = A*m\n\ntaper for a space containing 1D arrays, and only taper at the end\n\nA = JopTaper(JetSpace(Float64,10), (1,), (0.0,), (0.25,))\nm = ones(domain(A))\nd = A*m\n\ntaper for a space containing 2D arrays, where both dimensions are tapered\n\nA = JopTaper(JetSpace(Float64,10,11), (1,2), (0.75,0.5))\nm = ones(domain(A))\nd = A*m\n\ntaper for a space containing 3D arrays, where two of the three dimensions are tapered\n\nA = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.75,0.5))\nm = ones(domain(A))\nd = A*m\n\nIn the next example, array dimension 1 is tapered at the end, and array dimension 2 is tapered at the beginning and end.\n\nA = JopTaper(JetSpace(Float64,10,11,12), (1,3), (0.0,0.25), (0.25,0.25))\nm = ones(domain(A))\nd = A*m\n\n\n\n\n\n","category":"method"},{"location":"reference/#JetPack.JopTranslation-Union{Tuple{T}, Tuple{Matrix{T}, Matrix{T}}} where T","page":"Reference","title":"JetPack.JopTranslation","text":"A = JopTranslation(p1::Matrix, p2::Matrix)\n\nA*m translates m::Matrix using the vector field (p1,p2).\n\n\n\n\n\n","category":"method"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"#JetPack","page":"JetPack","title":"JetPack","text":"","category":"section"},{"location":"","page":"JetPack","title":"JetPack","text":"A somewhat hap-hazard set of operators for Jets.jl.","category":"page"},{"location":"#Example","page":"JetPack","title":"Example","text":"","category":"section"},{"location":"","page":"JetPack","title":"JetPack","text":"using Jets, JetPack, Random\nR = JetSpace(Float64,128,128)\nx = ones(R)\nA = JopRestriction(R, randperm(R, 512))\ny = A*x\nz = A'*y","category":"page"},{"location":"","page":"JetPack","title":"JetPack","text":"One can now plot x and z to see the effect of the restriction.  For example, using PyPlot:","category":"page"},{"location":"","page":"JetPack","title":"JetPack","text":"using PyPlot\nfigure(1);clf()\nsubplot(121);imshow(x)\nsubplot(122);imshow(z) ","category":"page"}]
}
