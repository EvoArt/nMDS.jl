module nMDS
using  Distances,  LoopVectorization,  DataStructures

include("functions.jl",
    "main.jl")
export nmds
end
