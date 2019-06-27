function sor(A, b;
    tol=sqrt(eps(real(eltype(b)))), # tolerance
    maxiter=length(b),              # maximum number of iterations
    omega=1.5
    )
    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(maxiter+1)
    res[1] = nb
    for k = 1:maxiter
        for i = 1:n
            sigma = 0.
            for j = 1:n
                if j != i
                    sigma += A[i,j] * x[j]
                end
            end
            sigma = (b[i] - sigma) / A[i,i]
            x[i] += omega * (sigma - x[i])
        end
        # Check residual and error
        res[k+1] = norm(A*x - b)
        print_residual(k, maxiter, res[k+1], nb)
        if res[k+1]/nb < tol
            success_message(k)
            return (x, res[1:k+1])
        end
    end
    failure_message(maxiter)
    return (x, res)
end
