using LinearAlgebra
using Printf

function gmres_update(x, s, q, i, H)
    y = H[1:i,1:i] \ s[1:i]
    for k in eachindex(y)
        x += q[k] * y[k]
    end
    return x
end

function gmres(A, b;
    tol = sqrt(eps(real(eltype(b)))), # tolerance
    restart = length(b), # maximum number of inner iterations between restarts
    maxiter = min(20, length(b)),     # maximum number of outer iterations
    precond = UniformScaling(1.0))    # preconditioner

    # It is assumed that the preconditioner supports the operator \
    # If you pass a matrix, we advise that you call lufact() before
    # calling this routine such that \ can be computed efficiently,
    # without having to refactor the matrix for every call of \.

    n = length(b)
    m = restart # variable with shorter name
    x = zeros(n)
    nb = norm(b)
    res = zeros(1 + maxiter * m)
    res[1] = nb

    M = precond # shorter name

    it  = 0

    for j = 1:maxiter

        residual = b - A * x
        q = Vector{Any}(undef, m + 1)
        J = Vector{Any}(undef, m)
        H = zeros(m + 1, m)

        @printf "Outer iteration %d/%d, residual %1.2e\n" j maxiter norm(residual)

        r = M \ residual
        normr = norm(r)
        q[1] = r / normr
        s = zeros(m + 1)
        s[1] = normr

        for i = 1:m
            z = M \ (A * q[i])
            # Arnoldi iteration
            for k = 1:i
                H[k,i] = z' * q[k]
                z -= H[k,i] * q[k]
            end
            H[i + 1,i] = norm(z)
            q[i + 1] = z / H[i + 1,i]

            # Apply previous Givens rotations to solve least squares
            for k = 1:i - 1
                H[1:i + 1,i] = J[k] * H[1:i + 1,i]
            end
            J[i], = givens(H[i,i], H[i + 1,i], i, i + 1)
            # Update s and H
            H[1:i + 1,i] = J[i] * H[1:i + 1,i]
            s = J[i] * s

            it += 1
            # Norm of residual
            res[it + 1] = abs(s[i + 1])
            print_residual(i, maxiter, res[it + 1], nb)

            # Check residual, compute x, and stop if possible
            if res[it + 1] / nb < tol
                x = gmres_update(x, s, q, i, H)
                @printf "Converged after %d outer and %d inner iterations\n" j i
                return (x, res[1:it + 1])
            end
        end

        # Update x before the restart
        x = gmres_update(x, s, q, m, H)
    end
    failure_message(maxiter)
    return (x, res)
end
