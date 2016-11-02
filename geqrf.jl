"""
    house(x)

Computes the Householder transformation for input vector x.
"""
function house(x)
    """Computes the Householder transformation for input vector x"""
    sigma = dot(x[2:end],x[2:end])
    v = copy(x)

    if sigma == 0
        beta = 0
        return beta, v
    end

    sq = sqrt(x[1]^2 + sigma)
    if x[1] > 0
        v[1] += sq
    else
        v[1] -= sq
    end

    beta = 2.0 / (v[1]^2 + sigma)

    return beta, v
end

"""
    geqrf!(A)

QR factorization of A using Householder transformation
"""
function geqrf!(A)
    m = size(A,1)
    n = size(A,2)
    vA = zeros(n)
    kend = (m > n ? n : m-1)
    for k=1:kend
        beta, v = house(A[k:end,k])
        for j=k:n
            # vA = beta * v^T * A
            vA[j] = 0.0
            for i=k:m
                vA[j] += v[i-k+1] * A[i,j]
            end
            vA[j] *= beta
        end
        # A - beta v (v^T A)
        for j=k:n, i=k:m
            A[i,j] -= v[i-k+1] * vA[j]
        end
        A[k+1:end,k] = v[2:end] 
        # Saving v in the lower triangular part of A.
        # This was not done here but one can always
        # divide v by v[1] such that v[1] = 1 is always true.
        # In that case, v[1] does not need to be stored.
    end
    # Lower triangular part of A: sequence of v vectors
    # Upper triangular part: factor R
end

"""
    givens(a, b)

Computes the Givens transform; a and b should be scalars
"""
function givens(a, b)
    if b == 0
        c = 1
        s = 0
    else
        if abs(b) > abs(a)
            tau = -a/b
            s = 1.0/sqrt(1.0+tau*tau)
            c = s*tau
        else
            tau = -b/a
            c = 1.0/sqrt(1.0+tau*tau)
            s = c*tau
        end
    end
    return (c, s)
end

"""
    geqrfCGS!(A)

QR factorization of A using Gram-Schmidt
"""
function geqrfCGS!(A)
    m = size(A,1)
    n = size(A,2)
    @assert m >= n
    R = zeros(Float64, n,n)
    for j=1:n
        # Orthogonalize
        for i=1:j-1, k=1:m
            R[i,j] += A[k,i] * A[k,j]
        end
        for i=1:j-1, k=1:m
            A[k,j] -= A[k,i] * R[i,j]
        end
        R[j,j] = norm( A[:,j] )
        # Normalize column
        A[:,j] /= R[j,j]
    end
    # A contains the Q factor at the end
    return R
end

"""
    geqrfMGS!(A)

QR factorization of A using modified Gram-Schmidt
"""
function geqrfMGS!(A)
    m = size(A,1)
    n = size(A,2)
    @assert m >= n
    R = zeros(Float64, n,n)
    for j=1:n
        # Orthogonalize
        for i=1:j-1
            for k=1:m
                R[i,j] += A[k,i] * A[k,j]
            end
            for k=1:m
                A[k,j] -= A[k,i] * R[i,j]
            end
        end
        # Normalize column
        R[j,j] = norm( A[:,j] )
        A[:,j] /= R[j,j]
    end
    # A contains the Q factor at the end
    return R
end
