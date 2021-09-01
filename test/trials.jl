using  Distances,  LoopVectorization, BenchmarkTools, DataStructures, GLMakie

# https://link.springer.com/content/pdf/10.1007/BF02289694.pdf
function stress(d,d̂,M)
    S_star = 0
    T_star = 0
    @turbo for m in 1:M
        S_star += (d[m] - d̂[m])^2
        T_star += d[m]^2
    end
    return S_star, T_star, sqrt(S_star/T_star)
end

function get_dist!(X1,DIST,Is,Js,M)
    @inbounds for m in 1:M
        DIST[m] = euclidean(X1[Is[m],:],X1[Js[m],:])
    end
end


# calculate the gratient g for K coordinates over L dimensions. 
# Since k and l each only appera in one term of the rhs, we can calculate K versions of
# the first term, one of the secod and L of the third. Then We can combine them to get g 

# the second (shared term)
function grhs2(d, d̂, S_star,T_star,M,res) 
    @turbo for m in 1:M
        res[m] = (d[m]-d̂[m])/S_star - d[m]/T_star
    end
    return res
end 
function grhs2!(d, d̂, S_star,T_star,M,res) 
    @turbo for m in 1:M
        res[m] = (d[m]-d̂[m])/S_star - d[m]/T_star
    end
end 


function grhs3!(x,d, M,L,res,Is,Js) 
    Imat = view(x,Is,:)
    Jmat = view(x,Js,:)
    @turbo for m in 1:M
        for l in 1:L
            res[m,l] = (Imat[m,l] - Jmat[m,l])/d[m]
        end
    end
end 

function g(grhs2,grhs3,L,M,S,Is,Js,X2,grhs)
    @inbounds for l in 1:L
       grhs .= grhs2 .* view(grhs3,:,l)
       for m in 1:M
           X2[Is[m],l] += grhs[m]
           X2[Js[m],l] -= grhs[m]
       end
   end
   return X2 *S
end




function monotone(d,M)
    loc = 1
    partition = collect(1:M)
    block = view(d,partition .==loc)
    #start at position 1, up-active i.e. check to see that the block to the right has a higher d̂
    while (loc < maximum(partition)) | (mean(block) < mean(view(d,partition .==loc-1)))
       # println(loc)
        block = view(d,partition .==loc)
    if (mean(block)< mean(view(d,partition .==loc+1))) | (loc == maximum(partition))
        if (mean(block)> mean(view(d,partition .==loc-1))) | (loc ==1)
            loc += 1
        else
            partition[partition .==loc] .=loc-1
            partition[partition .>=loc] .-=1
            loc -=1
        end
        else
            partition[partition .==loc+1] .=loc
            partition[partition .>loc] .-=1
        end
    end
    blockmax = partition[end]
    blocks = 1:blockmax
    
    d̂ = Vector{Float64}(undef,M)
    for block in blocks
        boolmask = partition .== block
        d̂[boolmask] .= mean(view(d,boolmask))
    end
    return d̂
end

function monotone!(d,d̂,M)
    loc = 1
    partition = collect(1:M)
    block = view(d,partition .==loc)
    #start at position 1, up-active i.e. check to see that the block to the right has a higher d̂
    while (loc < maximum(partition)) | (mean(block) < mean(view(d,partition .==loc-1)))
       # println(loc)
        block = view(d,partition .==loc)
    if (mean(block)< mean(view(d,partition .==loc+1))) | (loc == maximum(partition))
        if (mean(block)> mean(view(d,partition .==loc-1))) | (loc ==1)
            loc += 1
        else
            partition[partition .==loc] .=loc-1
            partition[partition .>=loc] .-=1
            loc -=1
        end
        else
            partition[partition .==loc+1] .=loc
            partition[partition .>loc] .-=1
        end
    end
    blockmax = partition[end]
    blocks = 1:blockmax
    
    for block in blocks
        boolmask = partition .== block
        d̂[boolmask] .= mean(view(d,boolmask))
    end
end


function normalize(x)
    x .-= mean(x)
    return x ./std(x)
end
function normalize(x,dim)
    x .-= mean(x,dims = dim)
    return x ./std(x,dims = dim)
end

cosθ(g,g″) = sum(g .* g″)/(sqrt(sum(g .^2))*sqrt(sum(g″ .^2)))
angle_factor(g,g″) = 4.0^(cosθ(g,g″)^3.0)
relaxation_factor(x) = 1.3/(1+(five_step_ratio(x)^5))
five_step_ratio(x) = min(1,x[end]/x[1])
good_luck_factor(x) = min(1,x[end]/x[end-1])
mag(g,K) = sqrt(sum(g .^2)/K)
x′(x,g,α,K) = x .+ g*α/mag(g,K)


K = 20
L = 2
X = rand(K,12)
X[1:10,:] .+= 0.5
D = pairwise(Euclidean(),X')
X1 = normalize(randn(K,L),1)
X2 = zeros(K,L)
X3 = zeros(K,L)
α = 0.2

stress_vec = CircularBuffer(6)
DISSIM = Vector{Float64}(undef,0)
IJ = Vector{Tuple}(undef,0)
Is = Vector{Int64}(undef,0)
Js = Vector{Int64}(undef,0)

for i in 2:K
    for j in 1:i-1 
        push!(DISSIM,D[i,j])
        push!(IJ,(i,j))
        push!(Is,i)
        push!(Js,j)
    end
end
M = length(IJ)
order = sortperm(DISSIM)
DISSIM = DISSIM[order]
IJ = IJ[order]
Is = Is[order]
Js = Js[order]
DIST = Vector{Float64}(undef,M)
for m in 1:M
    i,j = IJ[m]
    DIST[m] = euclidean(X1[i,:],X1[j,:])
end

DHAT = monotone(DIST,M)
S_star, T_star, S = stress(DIST,DHAT,M)
push!(stress_vec,S)
g2 = Vector{Float64}(undef,M)
g2=grhs2(DIST, DHAT, S_star,T_star,M,g2) 
g3 = Array{Float64}(undef,M,L)
g3 =grhs3(X1,DIST, M,L,g3,Is,Js) 
res = rand(K)
grhs = Vector{Float64}(undef,M)
X2 = g(g2,g3,L,M,S,Is,Js,X2,grhs) .*-1

X1 = normalize(x′(X1,X2,α,K),1)
X3 .= X2
X2 = zeros(K,L)

slist = []
for step in 1:150
for m in 1:M
    i,j = IJ[m]
    DIST[m] = euclidean(X1[i,:],X1[j,:])
end

DHAT = monotone(DIST,M)
S_star, T_star, S = stress(DIST,DHAT,M)
push!(stress_vec,S)
push!(slist,S)
g2=grhs2(DIST, DHAT, S_star,T_star,M,g2) 
g3 =grhs3(X1,DIST, M,L,g3,Is,Js) 

X2 = g(g2,g3,L,M,S,Is,Js,X2,grhs) .*-1

α = α * angle_factor(X2,X3) * relaxation_factor(stress_vec) * good_luck_factor(stress_vec)
X3 = X2
X1 = normalize(x′(X1,X2,α,K),1)
#X1 =x′(X1,X2,α,K)
X2 = zeros(K,L)
end
lines(float.(slist))
scatter(X1[1:10,1],X1[1:10,2])
scatter!(X1[11:20,1],X1[11:20,2], color = :red)


function nMDS(D,L,local_minimum_criterion = 0.02, max_iter = 10000)
    # initialise varables and containers
    local_minimum = false
    K = size(D)[1]
    X2 = zeros(K,L)
    X3 = zeros(K,L)
    α = 0.2
    rand_S = 0.0

    stress_vec = CircularBuffer(6)
    DISSIM = Vector{Float64}(undef,0)
    IJ = Vector{Tuple}(undef,0)
    Is = Vector{Int64}(undef,0)
    Js = Vector{Int64}(undef,0)
    @inbounds for i in 2:K
        for j in 1:i-1 
            push!(DISSIM,D[i,j])
            push!(IJ,(i,j))
            push!(Is,i)
            push!(Js,j)
        end
    end
    M = length(IJ)
    DIST = Vector{Float64}(undef,M)
    DHAT = Vector{Float64}(undef,M)
    order = sortperm(DISSIM)
    DISSIM = DISSIM[order]
    IJ = IJ[order]
    Is = Is[order]
    Js = Js[order]
    g2 = Vector{Float64}(undef,M)
    g3 = Array{Float64}(undef,M,L)
    grhs = Vector{Float64}(undef,M)
    X1 = normalize(randn(K,L),1)
    iter = 0
    while (local_minimum == false) & (iter <max_iter) 
        iter += 1
        get_dist!(X1,DIST,Is,Js,M)
        monotone!(DIST,DHAT,M)
        S_star, T_star, S = stress(DIST,DHAT,M)
        push!(stress_vec,S)   
        
        grhs2!(DIST, DHAT, S_star,T_star,M,g2)    
        grhs3!(X1,DIST, M,L,g3,Is,Js) 
        X2 = g(g2,g3,L,M,S,Is,Js,X2,grhs) .*-1
        α =  length(stress_vec) < 2 ? α : α * angle_factor(X2,X3) * relaxation_factor(stress_vec) * good_luck_factor(stress_vec)
        X1 = normalize(x′(X1,X2,α,K),1)
        X3 .= X2
        
        if length(stress_vec) == 1
            rand_S = mag(X2,K)
        elseif  mag(X2,K) < local_minimum_criterion*rand_S
            local_minimum = true
        end
    end
    return X1
end

   

