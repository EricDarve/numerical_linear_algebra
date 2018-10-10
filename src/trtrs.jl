function trtrs(L, b)
    n = length(b)
    x = copy(b)
    for j = 1:n
        x[j] = x[j] / L[j,j]
        for i = j+1:n
            x[i] -= L[i,j] * x[j]
            # This is fast because the row index is the inner loop
        end
    end
    return x
end

function trtrsRow(L, b) # This is the slow implementation
    n = length(b)
    x = Vector{Float64}(undef,n)
    for i=1:n
        z = 0.0
        for j = 1:i-1
            z += L[i,j] * x[j]
        end
        x[i] = (b[i] - z) / L[i,i]
    end
    return x
end
