#This function implements the algorithm from  https://link.springer.com/content/pdf/10.1007/BF02289694.pdf
# Varibales are generally named simillaryly to thos in the original paper. However
# it felt simples to have a a vector of Is and one of Js, rather than an IJ array.
# All functions used can be found in functions.jl

# The function takes a distance matrix D and returns a K,L array of transformed coordinates
# The algorithm exits once the magnitude of the gradient is less than the maginitude
# from a random arrangement of coordinates multiplied by `local_minimum_criterion`
# or once it has gone through `max_iter` iterations.

function nmds(D,L ::Int64,local_minimum_criterion = 0.02, max_iter ::Int64 = 10000)
    # initialise varables and containers
    local_minimum = false
    K = size(D)[1] # Number of observations/points
    X2 = zeros(K,L) # Current gradient
    X3 = zeros(K,L) # Previous gradient
    α = 0.2 # Initial step size
    rand_mag = 0.0 # To be updated as random magnitude of gradient
    stress_vec = CircularBuffer(6) # Current and previous 5 stress calculations
    DISSIM = Vector{Float64}(undef,0) # Original dissimilarities
    Is = Vector{Int64}(undef,0) # i indeces
    Js = Vector{Int64}(undef,0) # j indeces
    @inbounds for i in 2:K
        for j in 1:i-1 
            push!(DISSIM,D[i,j])
            push!(Is,i)
            push!(Js,j)
        end
    end
    M = length(Is) 
    DIST = Vector{Float64}(undef,M) # euclidean distances between points in X1
    DHAT = Vector{Float64}(undef,M) # fitted distances
    # gradient calculations will be broken into blocks to give vectors of length M
    g2 = Vector{Float64}(undef,M) # solutions to 2nd term on rhs of gradient calculation
    g3 = Array{Float64}(undef,M,L) # solutions to 3rd term on rhs of gradient calculation
    grhs = Vector{Float64}(undef,M) # combined 2nd and third terms for corresponding value of l
                                    # this will be repeatedly overwritten as we iterate over l
    # sort parameters by order of values in DISSIM
    order = sortperm(DISSIM)
    DISSIM = DISSIM[order]
    Is = Is[order]
    Js = Js[order]
    # generate random starting values for X1
    X1 = normalize(randn(K,L),1) # current point values
    iter = 0
    while (local_minimum == false) & (iter <max_iter) 
        iter += 1
        get_dist!(X1,DIST,Is,Js,M) # calculate euclidean distances between points
        monotone!(DIST,DHAT,M) # fit DHAT values vie monotone regression
        S_star, T_star, S = stress(DIST,DHAT,M) # calculate stress
        push!(stress_vec,S)   # accumulate stress in vector
        grhs2!(DIST, DHAT, S_star,T_star,M,g2) # calculate second term of gradient rhs
        grhs3!(X1,DIST, M,L,g3,Is,Js) # calculate 3rd ter of gradient rhs
        X2 = g(g2,g3,L,M,S,Is,Js,X2,grhs) # put them together and calculate negative gradient
        # if there are 2 or more stress calculations so far, use Kruskals method to get new step size
        α =  length(stress_vec) < 2 ? α : α * angle_factor(X2,X3) * relaxation_factor(stress_vec) * good_luck_factor(stress_vec)
        X1 = normalize(x′(X1,X2,α,K),1) # get new point values
        X3 .= X2 # update previous gradient
        
        if iter == 1 # on first iteration
            rand_mag = mag(X2,K) # magnitude of gradient
        elseif  mag(X2,K) < local_minimum_criterion*rand_mag # if local minimum 
            local_minimum = true # exit loop
        end
    end
    return X1#, S
end

# run the above algorith n times and return the result with the lowest stress value
function nmds(D,L ::Int64,n ::Int64, local_minimum_criterion = 0.02, max_iter ::Int64 = 10000)
    results = []
    Threads.@threads for i in 1:n
        push!(results,nmds(D,L,local_minimum_criterion,max_iter))
    end
    best = findmin([results[i][2] for i in 1:n])
    return results[best]
end
