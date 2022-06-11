function isotonic_regression!(y::Vector{Float64})
    n = length(y)
    @inbounds begin
        n -= 1
        while true
            i = 1
            pooled = 0
            while i <= n
                k = i
                while k <= n && y[k] >= y[k+1]
                    k += 1
                end

                # Find a decreasing subsequence, and update
                # all points in the sequence to the weighted average.
                if y[i] != y[k]
                    numerator = 0.0
                    #denominator = 0.0
                    @turbo for j in i : k
                        numerator += y[j] #* weights[j]
                        #denominator += 1 weights[j]
                    end
                    frac = numerator/((k - i) +1)
                    @turbo for j in i : k
                        y[j] = frac
                    end
                    pooled = 1
                end
                i = k + 1
            end
            if pooled == 0
                break
            end
        end
    end
    #return y
end

#isotonic_regression!(y::Vector{Float64}) = isotonic_regression!(y)

# non-mutating versions
isotonic_regression(y::Vector{Float64}) = isotonic_regression!(copy(y))
#isotonic_regression(y::Vector{Float64}) = isotonic_regression(y)

