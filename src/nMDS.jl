module nMDS
using  Distances,  LoopVectorization,  DataStructures

include("functions.jl")
include("main.jl")
export nmds
end
