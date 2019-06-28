using LinearAlgebra

function qmr(A, b;
    tol = sqrt(eps(real(eltype(b)))), # tolerance
    maxiter = length(b),              # maximum number of iterations
    precond1 = UniformScaling(1),     # preconditioner M1
    precond2 = UniformScaling(1))      # preconditioner M2

    # Preconditioner is defined as M = M1 M2.
    # It is assumed that the preconditioners M1 and M2 support the operator \
    # If you pass a matrix, we advise that you call lufact() before
    # calling this routine such that \ can be computed efficiently,
    # without having to refactor the matrix for every call of \.

    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(1 + maxiter)
    res[1] = nb

    M1 = precond1 # shorter name
    M2 = precond2 # shorter name

    vtilde = b
    y = M1 \ vtilde
    rho = norm(y)

    wtilde = b
    z = (M2') \ wtilde
    xi = norm(z)

    gamma = 1.
    eta = -1.

    r = b

    epsilon = 1.
    theta = 0.
    p = zeros(n)
    q = zeros(n)
    d = zeros(n)
    s = zeros(n)

    for i = 1:maxiter
        if rho == 0. || xi == 0.
            failure_message(i)
            return (x, res[1:i])
        end

        y /= rho
        z /= xi
        delta = dot(z, y)

        y = M2 \ y
        z = (M1') \ z
        p = y - (xi  * delta / epsilon) * p
        q = z - (rho * delta / epsilon) * q

        ptilde = A * p
        epsilon = dot(q, ptilde)
        if delta == 0. || epsilon == 0.
            failure_message(i)
            return (x, res[1:i])
        end
        beta = epsilon / delta

        vtilde = ptilde - (beta / rho) * vtilde
        y = M1 \ vtilde
        rhoold, rho = rho, norm(y)

        wtilde = A' * q - (beta / xi) * wtilde
        z = (M2') \ wtilde
        xi = norm(z)

        thetaold, theta = theta, rho / (gamma * abs(beta))
        gammaold, gamma = gamma, 1. / sqrt(1. + theta^2)

        eta *= - (rhoold / beta) / (gammaold / gamma)^2

        d = eta * p      + (thetaold * gamma)^2 * d
        s = eta * ptilde + (thetaold * gamma)^2 * s
        x += d
        r -= s

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
