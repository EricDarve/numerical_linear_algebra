using LinearAlgebra

function cr(A, b;
    tol = sqrt(eps(real(eltype(b)))), # tolerance
    maxiter = length(b),              # maximum number of iterations
    precond = UniformScaling(1.0))    # preconditioner

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

    r = M \ b

    rho = 1.
    p = zeros(n)
    q = p

    for i = 1:maxiter
        s = A * r
        rhoold, rho = rho, dot(r, s)
        beta = rho / rhoold
        p = r + beta * p
        q = s + beta * q
        z = M \ q
        alpha = rho / dot(q, z)
        x += alpha * p
        r -= alpha * z
        res[i + 1] = norm(r)
        print_residual(i, maxiter, res[i + 1], nb)
        if res[i + 1] / nb < tol
            success_message(i)
            return (x, res[1:i + 1])
        end
    end
    failure_message(maxiter)
    return (x, res)
end
