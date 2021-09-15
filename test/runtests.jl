using Revise, ProfileView
using nMDS, Distances
using Test
x = rand(10,10)
dist = pairwise(BrayCurtis(),x)
@profview nmds(dist,3)
@profview nmds(dist,3)



@testset "nMDS.jl" begin
    # Write your tests here.
end
