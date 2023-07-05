module JetPackFluxExt

using JetPack, Jets, Flux
"""
    F = JopAD(f, dom, rng)

where `F` wraps a Julia native function `f` as a Jet operator with domain
and range given by `dom` and `rng`, respectively. Notice that it is the
responsibility of the user to ensure that function `f` is auto-differentiable
in both forward-mode and reverse-mode via the function `pushforward` and
`pullback` in Flux. If the user chooses to use customized AD rules for this
function `f`, it is also the responsibility of the user to ensure that these
rules are tested (e.g. linearity, dot product, and linearization tests).
""" 
JetPack.JopAD(f::Function; dom, rng) = JopNl(dom = dom, rng = rng, f! = JopAD_f!, df! = JopAD_df!, df′! = JopAD_df′!, s = (f=f,))

JopAD_f!(d::AbstractArray, m::AbstractArray; f::Function) = d .= f(m)
JopAD_df!(δd::AbstractArray, δm::AbstractArray; mₒ, f::Function) = δd .= Flux.pushforward(f, mₒ)(δm)
JopAD_df′!(δm::AbstractArray, δd::AbstractArray; mₒ, f::Function) = δm .= Flux.pullback(f, mₒ)[2](δd)[1]

end