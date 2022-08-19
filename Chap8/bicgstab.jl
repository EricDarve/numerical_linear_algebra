using LinearAlgebra

function bicgstab(A, b;
    tol = sqrt(eps(real(eltype(b)))), # tolerance
    maxiter = length(b),              # maximum number of iterations
    precond = UniformScaling(1.0))     # preconditioner

    # It is assumed that the preconditioner supports the operator \
    # If you pass a matrix, we advise that you call lufact() before
    # calling this routine such that \ can be computed efficiently,
    # without having to refactor the matrix for every call of \.

    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(maxiter + 1)
    res[1] = nb

    r = b
    rtilde = b

    M = precond # shorter name

    rho = 1.
    alpha = 0.
    omega = 1.
    p = zeros(n)
    v = p

    for i = 1:maxiter
        rhoold, rho = rho, dot(rtilde, r)
        if rho == 0. || omega == 0.
            failure_message(i)
            return (x, res[1:i])
        end
        p = r + (rho / rhoold) * (alpha / omega) * (p - omega * v)
        phat = M \ p
        v = A * phat
        alpha = rho / dot(rtilde, v)
        s = r - alpha * v
        if norm(s) / nb < tol
            res[i + 1] = norm(s)
            print_residual(i, maxiter, res[i + 1], nb)
            success_message(i)
            x += alpha * phat
            return (x, res[1:i + 1])
        end
        shat = M \ s
        t = A * shat
        omega = dot(t, s) / dot(t, t)
        x += alpha * phat + omega * shat
        r = s - omega * t
        res[i + 1] = norm(r)
        @show res[i + 1] / nb, tol
        print_residual(i, maxiter, res[i + 1], nb)
        if res[i + 1] / nb < tol
            success_message(i)
            return (x, res[1:i + 1])
        end
    end
    failure_message(maxiter)
    return (x, res)
end
