using LinearAlgebra

function symortho(a, b)
    if (b == 0)
        s = 0
        r = abs(a)
        if a == 0
            c = 1
        else
            c = sign(a)
        end
    elseif a == 0
        c = 0
        r = abs(b)        
        s = sign(b)
    elseif abs(b) > abs(a)
        tau = a / b
        s = sign(b) / sqrt(1 + tau^2)
        c = s * tau
        r = b / s
    else
        tau = b / a
        c = sign(a) / sqrt(1 + tau^2)
        s = c * tau
        r = a / c
    end
    return (c, s, r)
end

function minres(A, b;
    tol = sqrt(eps(real(eltype(b)))), # tolerance
    maxiter = length(b),              # maximum number of iterations
    precond = UniformScaling(1.0))     # preconditioner

    println(".............V4.............")

    # It is assumed that the preconditioner supports the operator \
    # If you pass a matrix, we advise that you call lufact() before
    # calling this routine such that \ can be computed efficiently,
    # without having to refactor the matrix for every call of \.

    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(maxiter + 1)
    res[1] = nb

    M = precond # shorter name

    q = M \ b
    beta = norm(q)
    q /= beta
    qold = zeros(n)

    phi = beta

    c, s = -1., 0.
    delta = 0.
    
    eps = 0.
    pold, p = zeros(n), zeros(n)

    for i = 1:maxiter
        # Lanczos step
        z = M \ (A * q)
        alpha = dot(q, z) # Diagonal entry in Tk
        qold, q = q, z - alpha * q - beta * qold

        beta = norm(q) # Off-diagonal entry in Tk
        if beta == 0.
            failure_message(i)
            return (x, res[1:i])
        end

        q /= beta # Normalize q vector

        # Solve min ||Tk' y - beta0 e1||
        tau  = c * delta + s * alpha
        gamma = s * delta - c * alpha

        delta = - c * beta
        epsold, eps = eps, s * beta

        # Build rotation
        c, s, eta = symortho(gamma, beta)

        mu = c * phi
        phi = s * phi

        # Update search direction p
        pold, p = p, (qold - tau * p - epsold * pold) / eta

        # Update solution
        x += mu * p

        # Check residual
        res[i + 1] = abs(phi)
        print_residual(i, maxiter, res[i + 1], nb)
        if res[i + 1] / nb < tol
            success_message(i)
            return (x, res[1:i + 1])
        end
    end
    failure_message(maxiter)
    return (x, res)
end
