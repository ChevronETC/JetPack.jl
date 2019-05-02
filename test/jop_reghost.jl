using Revise
using JetPack, Jets, Test

@testset "reghost" begin
    nt = 512
    nx = 256
    G = JopReghost(symspace(JopReghost, Float64, nt, nx), 10.0, 20.0, .004, 10.0) # target depth, receiver depth, time sampling, space sampling
    m = rand(domain(G))
    d = rand(range(G))
    lhs, rhs = dot_product_test(G, m, d)
    @test lhs ≈ rhs
end
