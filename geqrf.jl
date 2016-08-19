function house(x)
    """Computes the Householder transformation for input vector x"""
    sigma = dot(x[2:end],x[2:end])
    v = copy(x)
    v[1] = 1.0
    if sigma == 0 && x[1] >= 0
        # x is already pointing in the +e_1 direction
        # We do nothing
        beta = 0.0
    elseif sigma == 0 && x[1] < 0
        # x is pointing in the -e_1 direction
        # This is a reflection with normal vector e_1
        beta = -2.0
    else
        # mu = ||x||
        mu = sqrt(x[1]*x[1] + sigma)
        if x[1] <= 0
            # Simple formula x_1 - ||x||
            v[1] = x[1] - mu
        else
            # Special formula for x[1] > 0 case
            v[1] = -sigma / (x[1] + mu)
        end
        # We set v[1] = 1 below so beta is given by
        # beta = 2 v[1]^2 / ||v||^2
        beta = 2.0 * v[1]^2 / (sigma + v[1]^2)
        # This sets v[1] = 1; this simplifies the storage of v
        # later on. Only entries v[2:end] will need to be stored.
        v /= v[1]
    end

    return beta, v
end

function geqrf(A)
    """QR factorization of A using Householder transformation"""
    n = size(A,2)
    for k=1:n-1
        v, beta = house(A[k:end,k])
        A[k:end,k:end] = A[k:end,k:end] - beta * dot(v, A[k:,k:]) * v
        A[k+1:end,k] = v # Saving v in the lower triangular part of A
    end
    # Lower triangular part of A: sequence of v vectors
    # Upper triangular part: factor R
    return A
end
