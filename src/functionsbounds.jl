# Helper functions for nMDS based on # https://link.springer.com/content/pdf/10.1007/BF02289694.pdf
# Varibales are generally named simillaryly to thos in the original paper. However
# it felt simples to have a a vector of Is and one of Js, rather than an IJ array.
# The complete algorithm, using these functions can be found in main.jl


#### Gradient calculations (for the Euclidean case) pp.126 of Kruskal(1964)

# There are 3 terms within the sum on the right hand side of the equation
# The first term simply denotes which values affect which elements of X2 and
# whether to add or subtract them. The second and third terms can be calculated
# separately to form a length M vector and an M by L arrau respectively. Then we
# loop over L and M, using the indeces in Is and Js to add t appropriate values to X2.
# We multiply by minus S, since we are interested in the negative gradient.

# Second term
function grhs2!(d, d̂, S_star,T_star,M,res) 
    @turbo for m in 1:M
        res[m] = (d[m]-d̂[m])/S_star - d[m]/T_star
    end
end 
# Third term
function grhs3!(x,d, M,L,res,Is,Js) 
    Imat = view(x,Is,:)
    Jmat = view(x,Js,:)
    @turbo for m in 1:M
        for l in 1:L
            res[m,l] = (Imat[m,l] - Jmat[m,l])/d[m]
        end
    end
end 
# put them together
function g(grhs2,grhs3,L,M,S,Is,Js,X2,grhs)
    @inbounds for l in 1:L
       grhs .= grhs2 .* view(grhs3,:,l)
       for m in 1:M
           X2[Is[m],l] += grhs[m]
           X2[Js[m],l] -= grhs[m]
       end
   end
   return X2 *-S
end

function checkup(start,block,partition,d,M)
    println("up")

    ∑ = 0
    i = start
    while i < M
    #while partition[i] == block
        println(i,"","hi")
        if partition[i] !== block
            break
        end
        ∑+=d[i]
        i+=1
    end
    return ∑/(i -start), i-1
end
function checkdown(start,block,partition,d)
    println("down")
    ∑ = 0
    i = start
    while  i>0
        println(i,"","hip")
        if partition[i] !== block
            break
        end
        ∑+=d[i]
        i-=1
        
    end
    return ∑/(start-i), i+1
end


#### Monotone rgression pp.126-128 Kruscal(1964)
#The R implementation monMDS also explains this procedure https://github.com/vegandevs/vegan/blob/master/src/monoMDS.f

# The aim is to sort the distances (in DIST) into blocks, such that the mean (d̂) of each block increases
# from first to last i.e. for any given block - the block before it will have a lower d̂ and the one after it will
# have a greater d̂. In Kruskals terms, each block must be downsatisfied (block before it has lower d̂) and up 
# satisfied (block after it has higher d̂).
# To achieve this we start with each element of d having it's own singleton block.
# Starting with the first block, we alternate between checking whether the block is upsatisfied and downsatisfied, whenever the block 
# is not satisfied for a given direction we merge it with the adjacent block in that direction and begin 
# checking that this block is up and down satified. If the current block is both up and down satisfied,
# we move on to the next block and do the same. This procedure continues until all blocks are up and down
# satisfied. 
# The fitted values returned to DHAT are the d̂ values associated with the block each DIST value belongs to. 

function monotone!(d,d̂,M, partition = false)
    loc = 1
    # if a starting partition is supplied e.g. because of ties, use it
    partition = partition == false ? collect(1:M) : partition
    block = view(d,partition .==loc)
    localmean = d[1]
    bounds = (1,1)
    #start at position 1, up-active i.e. check to see that the block to the right has a higher d̂
    while true#(loc < maximum(partition)) | (localmean < checkdown(bounds[1],loc-1,partition,d)[1])
    println(loc)
    println(bounds)
    println(partition[1:10])
    if loc < maximum(partition)   
        upmean, upbound = checkup(bounds[2]+1,loc+1,partition,d,M)
        println(localmean,"",upmean)
        if upmean > localmean # if up satissfied
            if loc > 1 # and not first block
                downmean, downbound = checkdown(bounds[1]-1,loc-1,partition,d)#check down
                if downmean < localmean # if down satisfied also
                    loc += 1 # move to next block
                    localmean = upmean
                    bounds = (bounds[2]+1,upbound) 

                else # if not down satisfied
                    partition[bounds[1]:end] .-=1 
                    localmean = (localmean*(bounds[2] -bounds[1]+1)+downmean*(bounds[1]-downbound))/(bounds[2]-downbound+1)
                    bounds = (downbound,bounds[2]) # merge down
                    loc -= 1
                end
            else #if first block and up satisfied
                loc += 1 # move to next block
                localmean = upmean
                bounds = (bounds[2]+1,upbound) 

            end
        else # if not upsatisfied
            partition[bounds[2]+1:end] .-=1 
            localmean = (localmean*(bounds[2] -bounds[1]+1)+upmean*(upbound-bounds[2]))/(upbound-bounds[1]+1)
            bounds = (bounds[1],upbound) # merge up
            #loc += 1
            if loc > 1 # and not first block
                downmean, downbound = checkdown(bounds[1]-1,loc-1,partition,d)#check down
              #  if downmean < localmean # if down satisfied also
                 #   loc += 1 # move to next block
                if   downmean > localmean  
               # else # if not down satisfied
                    partition[bounds[1]:end] .-=1 
                    localmean = (localmean*(bounds[2] -bounds[1]+1)+downmean*(bounds[1]-downbound))/(bounds[2]-downbound+1)
                    bounds = (downbound,bounds[2]) # merge down
                    loc -= 1
                 end
            end
        end
    else # if last block 
        println("last block")
        downmean, downbound = checkdown(bounds[1]-1,loc-1,partition,d)#check down
        if downmean < localmean # if down satisfied
            break # exit loop
        else # if not down satisfied
            partition[bounds[1]:end] .-=1 
            localmean = (localmean*(bounds[2] -bounds[1]+1)+downmean*(bounds[1]-downbound))/(bounds[2]-downbound+1)
            bounds = (downbound,bounds[2]) # merge down
            loc -= 1
        end
    end

    end

    blockmax = partition[end]
    blocks = 1:blockmax
    
    for block in blocks
        boolmask = partition .== block
        d̂[boolmask] .= mean(d[boolmask])
    end
end


# normalize point coordinates
function normalize(x,dim)
    x .-= mean(x,dims = dim)
    return x ./std(x,dims = dim)
end

#### Step size related function pp.120-122 of Kruskal(1964)
cosθ(g,g″) = sum(g .* g″)/(sqrt(sum(g .^2))*sqrt(sum(g″ .^2)))
angle_factor(g,g″) = 4.0^(cosθ(g,g″)^3.0)
relaxation_factor(x) = 1.3/(1+(five_step_ratio(x)^5))
five_step_ratio(x) = min(1,x[end]/x[1])
good_luck_factor(x) = min(1,x[end]/x[end-1])
mag(g,K) = sqrt(sum(g .^2)/K)
x′(x,g,α,K) = x .+ g*α/mag(g,K)

# Calculate the stress (S), S⋆ an T⋆ according to Kruskals stress formula pp.115 & 125 Kruskal(1964)
function stress(d,d̂,M)
    S_star = 0
    T_star = 0
    @turbo for m in 1:M
        S_star += (d[m] - d̂[m])^2
        T_star += d[m]^2
    end
    return S_star, T_star, sqrt(S_star/T_star)
end
# calculate Euclidean distances between points
function get_dist!(X1,DIST,Is,Js,M)
    @inbounds for m in 1:M
        DIST[m] = euclidean(X1[Is[m],:],X1[Js[m],:])
    end
end

