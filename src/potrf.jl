function potrf!(A)
    n = size(A,1)
    for j=1:n
        for k=1:j-1, i=j:n
            A[i,j] -= A[i,k] * A[j,k]
        end
        ajj = sqrt(A[j,j])
        for i=j:n
            A[i,j] /= ajj
        end
    end
end

function potrs(A,x)
    n = length(b)
    x = copy(b)
    # Solve using matrix G
    for j = 1:n
        x[j] /= A[j,j]
        for i = j+1:n
            x[i] -= A[i,j] * x[j]
        end
    end
    # Solve using matrix G^T
    for i = n:-1:1
        for j = i+1:n
            x[i] -= A[j,i] * x[j]
        end
        x[i] /= A[i,i]
    end
    return x
end
