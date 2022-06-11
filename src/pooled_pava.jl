function pooled_pava_isotonic_regression!(y::Vector{Float64},tievec = [], tiedict = Dict())

    n = length(y)
    j = 1
    S = Dict(0 => 0, 1 => 1)
    ydash = Dict(1 => y[1])
    wdash = Dict(1 => 1.0)

    i = 2
    @inbounds begin
        while i <= n
            iplus = 1
            j = j + 1
            #if in(i,tievec)
             #   tiedict_i = tiedict[i]
              #  iplus = tiedict_i-i+1
               # ydash[j] = mean(y[i:tiedict_i])
               # wdash[j] = sum(iplus)

            #else
                ydash[j] = y[i]
                wdash[j] = 1.0
           # end
            while j > 1 && ydash[j] < ydash[j-1]
                ydash[j-1] = (wdash[j] * ydash[j] + wdash[j-1] * ydash[j-1]) /
                   (wdash[j] + wdash[j-1])
                wdash[j-1] = wdash[j] + wdash[j-1]
                j = j-1
            end
            S[j] = i+iplus -1
            i += iplus
        end
        for k in 1 : j
            for l in S[k-1] + 1 : S[k]
                y[l] = ydash[k]
            end
        end
    end
    #return y
end

#pooled_pava_isotonic_regression!(y::Vector{Float64}) = pooled_pava_isotonic_regression!(y, ones(size(y, 1)))

# non-mutating versions
#pooled_pava_isotonic_regression(y::Vector{Float64}, weights::Vector{Float64}) = pooled_pava_isotonic_regression!(copy(y), weights)
pooled_pava_isotonic_regression(y::Vector{Float64}) = pooled_pava_isotonic_regression!(y, ones(size(y,1)))