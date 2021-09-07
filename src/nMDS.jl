module nMDS
using  Distances,  LoopVectorization,  DataStructures, Statistics

include("functions.jl")
include("main.jl")
export nmds
end
