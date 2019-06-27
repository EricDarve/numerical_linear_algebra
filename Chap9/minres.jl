function minres(A, b;
    tol=sqrt(eps(real(eltype(b)))), # tolerance
    maxiter=length(b),              # maximum number of iterations
    precond=speye(Float64, size(A,1), size(A,2)) # preconditioner
    )
    # It is assumed that the preconditioner supports the operator \
    # If you pass a matrix, we advise that you call lufact() before
    # calling this routine such that \ can be computed efficiently,
    # without having to refactor the matrix for every call of \.

    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(maxiter+1)
    res[1] = nb

    M = precond # shorter name

    v = M \ b
    beta = norm(v)
    v /= beta
    vold = zeros(n)

    phi = beta

    c, s = -1., 0.
    delta = 0.
    epsilon = 0.
    dold, d = 0., 0.

    for i = 1:maxiter
        # Lanczos step
        p = M \ (A*v)
        alpha = dot(v,p)
        p = p - alpha*v
        vold, v = v, p - beta*vold

        beta = norm(v)
        if beta == 0.
            failure_message(i)
            return (x, res[1:i])
        end

        v /= beta

        # Solve min ||Tk' z - beta0 e1||
        mu    = c * delta + s * alpha
        gamma = s * delta - c * alpha

        delta = - c * beta
        epsilonold, epsilon = epsilon, s * beta

        # Build next rotation
        eta = sqrt(gamma^2 + beta^2)
        c = gamma / eta
        s = beta  / eta

        tau = c * phi
        phi = s * phi

        # Update search direction d
        dold, d = d, (vold - mu * d - epsilonold * dold) / eta

        # Update solution
        x += tau * d

        # Check residual
        res[i+1] = abs(phi)
        print_residual(i, maxiter, res[i+1], nb)
        if res[i+1]/nb < tol
            success_message(i)
            return (x, res[1:i+1])
        end
    end
    failure_message(maxiter)
    return (x, res)
end
