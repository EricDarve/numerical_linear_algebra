function pretty_print(A)
    m = size(A,1)
    n = size(A,2)
    for i=1:m
        for j=1:n
            @printf "%9.1e" A[i,j]
        end
        @printf "\n"
    end
end

function lt(z1,z2)
    if abs( abs(z1)-abs(z2) ) > eps(Float32) * (abs(z1) + abs(z2))
        return abs(z1) < abs(z2)
    elseif abs( real(z1-z2) ) > eps(Float32) * (abs(z1) + abs(z2))
        return real(z1) < real(z2)
    else
        return imag(z1) < imag(z2)
    end
end

"""
    eigenvalue_sorted(A)

Sort eigenvalues of A first by the real part then by the imaginary part.
Assumes that A is real and we have pairs of complex conjugate eigenvalues.
"""
function eigenvalue_sorted(A)
    F = eigen(A)
    D = F.values
    sort!(D, lt=lt)
    return D
end

"""
    house(x)

Computes the Householder transformation for input vector x
"""
function house(x)
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
    apply_left_householder!(A, row, col, beta, v)

Apply a Householder reflection to the left of matrix `A`.

`row` is the starting row index to apply the transform.
`col` is the starting column index. `beta` and `v` are the
parameters for the Householder reflection.

"""
function apply_left_householder!(A, row, col, beta, v)
    n = size(A,1)
    vA = zeros(n)
    lv = length(v)
    # Apply transform to the left
    # vA = beta * v^T * A
    for j=col:n
        for i=row:row+lv-1
            vA[j] += v[i-row+1] * A[i,j]
        end
        vA[j] *= beta
    end
    # A - beta v (v^T A)
    for j=col:n, i=row:row+lv-1
        A[i,j] -= v[i-row+1] * vA[j]
    end
end

"""
    apply_right_householder!(A, row, col, beta, v)

Apply a Householder reflection to the right of matrix `A`.

`row` is the ending row index to apply the transform.
`col` is the starting column index. `beta` and `v` are the
parameters for the Householder reflection.
"""
function apply_right_householder!(A, row, col, beta, v)
    n = size(A,1)
    Av = zeros(n)
    lv = length(v)
    # Apply transform to the right
    # Av = beta * A * v
    for j=col:col+lv-1
        for i=1:row
            Av[i] += v[j-col+1] * A[i,j]
        end
    end
    for i=1:row
        Av[i] *= beta
    end
    # A - beta (Av) v^T
    for j=col:col+lv-1, i=1:row
        A[i,j] -= Av[i] * v[j-col+1]
    end
end

function gehrd!(A)
    """Reduced square matrix A to upper Hessenberg form"""
    n = size(A,1)
    for k=1:n-2
        # Compute Householder reflection for column k
        beta, v = house(A[k+1:n,k])
        # Apply the reflection on the left
        apply_left_householder!(A, k+1, k, beta, v)
        # Apply the reflection on the right
        apply_right_householder!(A, n, k+1, beta, v)
    end
end

function apply_left_givens!(A, row, col, c, s)
    n = size(A,1)
    for j=col:n
        A[row,j], A[row+1,j] =
            ( c * A[row,j] - s * A[row+1,j],
              s * A[row,j] + c * A[row+1,j] )
    end
end

function apply_right_givens!(A, row, col, c, s)
    n = size(A,1)
    for i=1:row
        A[i,col], A[i,col+1] =
            ( c * A[i,col] - s * A[i,col+1],
              s * A[i,col] + c * A[i,col+1] )
    end
end

"""
    givens_QR_iteration_s!(A)

QR iteration using Givens rotations.
"""
function givens_QR_iteration_s!(A)
    """First version of the algorithm to apply Givens for the QR iteration"""
    n = size(A,1)
    G = zeros(2,n-1)
    for k=1:n-1
        c, s = givens(A[k,k], A[k+1,k])
        G[:,k] = [c; s]
        # Multiply by Givens rotation to the left
        apply_left_givens!(A, k, k, c, s)
    end
    for k=1:n-1
        # Multiply by Givens rotation to the right
        apply_right_givens!(A, k+1, k, G[1,k], G[2,k])
    end
end

"""
    givens_QR_iteration!(A)

QR iteration with bulge chasing
"""
function givens_QR_iteration!(A)
    n = size(A,1)
    for k=1:n-1
        c, s = givens(A[k,k], A[k+1,k])
        apply_left_givens!(A, k, max(1,k-1), c, s)
        apply_right_givens!(A, min(n,k+2), k, c, s)
    end
end

"""
    double_shift_st(A)

For a real 2x2 matrix, calculate the sum and product of its eigenvalues.
"""
function double_shift_st(A)
    a = A[1,1]
    b = A[1,2]
    c = A[2,1]
    d = A[2,2]
    s = a+d       # sum
    t = a*d - b*c # product
    return s,t
end

"""
    gees_single_step!(A)

Single step in QR iteration with double real shift
"""
function gees_single_step!(A, exceptional_shift)
    n = size(A,1)

    @assert n>=3

    tol = 0.01
    # This tolerance is used to test for early convergence of the last
    # 2x2 or 1x1 block.

    # Which shift should we apply?
    if abs(A[n-1,n-2]) < tol * (abs(A[n-2,n-2]) + abs(A[n-1,n-1])) ||
   ! ( abs(A[n,n-1])   < tol * (abs(A[n-1,n-1]) + abs(A[n,n]))     )
        # Either:
        # (double shift test) == true
        # or
        # ! (single shift test) == true
        # double shift is the fall back position if the single shift test fails.
        s, t = double_shift_st(A[n-1:n,n-1:n])

    else # Single shift should be used
        s = 2*A[n,n]
        t = A[n,n]^2
    end

    if exceptional_shift
        println("Exceptional shift")
        # This case is difficult to handle.
        # A special shift is required to break the convergence
        # stall. We use the reference eig function for this.
        # Of course, this is not a "practical" solution.
        F = eigen(A); D = F.values
        D = sort(D, lt=lt)
        root = D[1]
        s = 2.0*real(root)
        t = abs(root)^2
    end

    # Assembling the first column
    v = [ A[1,1]*A[1,1] + A[1,2]*A[2,1] - s*A[1,1] + t;
          A[2,1]*(A[1,1]+A[2,2]-s);
          A[2,1]*A[3,2] ]

    beta, v = house(v)

    apply_left_householder!( A, 1, 1, beta, v)
    apply_right_householder!(A, min(4,n), 1, beta, v)

    for k=2:n-1
        beta, v = house(A[k:min(k+2,n),k-1])
        apply_left_householder!( A, k, k-1, beta, v)
        apply_right_householder!(A, min(k+3,n), k, beta, v)
    end
end

function reduce_eps!(A, tol)
    # Zero out all small entries on the sub-diagonal
    n = size(A,1)
    for i=1:n-1
        if abs(A[i+1,i]) < tol * (abs(A[i,i])+abs(A[i+1,i+1]))
            A[i+1,i] = 0
        end
    end
end

function gees!(A)
    n = size(A,1)
    if n==1
        return fill(A[1,1],1)
    end
    D = zeros(Complex{Float64},n)

    # Tolerance for deflation
    tol = eps(Float64)

    q = n # Size of the matrix we are currently working with
    iter = 1 # Counter to detect convergence failure
    iter_per_evalue = 0 # Used to trigger an exceptional shift

    # Zero out small entries
    reduce_eps!(A, tol)

    while q > 0

        if iter > 10*n
            println("Failure to converge")
            println("Size of matrix = ",n)
            println("Eigenvalues found:")
            println(D)
            D1 = eigenvalue_sorted(A)
            println("Remaining eigenvalues at point of failure:")
            println(D1)
            return D # The eigenvalues we were able to calculate up to the failure
        end

        deflation = true # Were we able to deflate the matrix?

        while deflation
            deflation = false
            # Test for convergence
            # Test deflation for the last 2x2 block
            if q <= 2 || A[q-1,q-2] == 0
                if q >= 2
                    # The last 2x2 block has converged
                    deflation = true # Deflating now
                    # Compute the eigenvalues
                    a = A[q-1,q-1]; b = A[q-1,q]; c = A[q,q-1]; d = A[q,q]
                    # Last 2x2 block
                    htr = (a+d)/2         # Half-trace
                    dis = (a-d)^2/4 + b*c # Discriminant
                    if dis > 0 # Pair of real eigenvalues
                        D[q-1] = htr - sqrt(dis)
                        D[q]   = htr + sqrt(dis)
                    else # Complex conjugate eigenvalues
                        D[q-1] = htr - sqrt(-dis)*im
                        D[q]   = htr + sqrt(-dis)*im
                    end
                    # Reduce the size of the matrix
                    q -= 2
                    if q>=1
                        A = A[1:q,1:q]
                    end
                end
            end

            if q==0
                return D
            end

            # Testing deflation for the last 1x1 block
            if q <= 1 || A[q,q-1] == 0
                deflation = true
                D[q] = A[q,q]
                q -= 1
                if q>=1
                    A = A[1:q,1:q]
                end
            end

            if q==0
                return D
            end

            if deflation
                iter_per_evalue = 0 # Reset the counter
            end
        end

        # If q <= 2 we will compute the eigenvalues at the next iteration
        if q >= 3
            # Searching for the smallest unreduced sub-block
            p = q
            while p > 1 && A[p,p-1] != 0
                p -= 1
            end

            # If the unreduced sub-block has size 2 or less, we move on
            # to the next iteration.
            if q-p+1 >= 3
                B = A[p:q,p:q] # Extract sub-block
                exceptional_shift = ((iter_per_evalue%5) == 0 &&
                                      iter_per_evalue>0)
                # Francis QR step
                gees_single_step!(B, exceptional_shift)
                # Reduce matrix
                reduce_eps!(B, tol)
                A[p:q,p:q] = B # Copy the resulting matrix back
                iter += 1 # Increment iteration counter
                iter_per_evalue += 1
                # Increment counter for exceptional_shift
            end
        end
    end
end

using Test
function gees_testsuite()
    n_case = 5       # Number of test cases
    s_test = 10      # Length of each test case

    err = zeros(n_case*s_test)
    h = 1e3*eps(Float64)

    itest = 1  # Test counter

    for i_case = 1:n_case
        for t = 1:s_test
            println("\n *** Test no ",itest," ***")

            rng = MersenneTwister(2016+itest)

            if i_case == 1 && t <= 4
                n = t
                A = rand(rng, n, n) # Small random
            elseif i_case == 1 && t == 5
                # Demmel case 1, p. 173
                A = Float64[0 0 1.0; 1.0 0 0; 0 1.0 0]
            elseif i_case == 1 && t == 6
                # Demmel case 2, p. 173
                A = Float64[0 1 0 0; 1 0 h 0; 0 -h 0 1; 0 0 1 0]
            elseif i_case == 1
                # Random matrix
                n = t
                A = zeros(n,n)
                n0 = div(n,2)
                n1 = n - n0
                As = rand(rng, n0, n0)
                A[1:n0,1:n0] = As + As'
                Ass = rand(rng, n1, n1)
                A[n-n1+1:n,n-n1+1:n] = Ass - Ass' + rand(rng) * UniformScaling(1.0)
                X = rand(rng, n, n);
                A = X * A / X
            elseif i_case == 2
                # Simple permutation
                n = t
                A = zeros(n,n)
                for j=1:n, i=1:n
                    if i == j+1
                        A[i,j] = 1
                    elseif i==1 && j==n
                        A[i,j] = 1
                    end
                end
            elseif i_case == 3
                # Random permutation
                # Eigenvalues are on the unit circle
                n = t
                r_perm = randperm(rng, n)
                A = zeros(n,n)
                for j=1:n, i=1:n
                    if r_perm[i] == j
                        A[i,j] = 1
                    end
                end
            elseif i_case == 4
                # Random permutation with small perturbation
                n = t
                r_perm = randperm(rng, n)
                A = zeros(n,n)
                for j=1:n, i=1:n
                    if r_perm[i] == j
                        A[i,j] = 1
                    end
                end
                if n >= 2
                    i = 1
                    if A[i,1] != 0
                        i = 2
                    end
                    A[i,1] = h
                end
            else
                # Random orthogonal matrix
                # Eigenvalues are on the unit circle
                n = t
                A = rand(rng, n, n)
                F = qr(A); A = Matrix(F.Q)
            end

            n = size(A,1)
            println("Size of matrix ",n)

            D0 = eigenvalue_sorted(A)

            # Reduction to upper Hessenberg form
            gehrd!(A)

            if n >= 3
                @test norm(tril(A,-2)) < 1e-8
                gees_single_step!(A,false)
                D1 = eigenvalue_sorted(A)
                @test norm(D0 - D1) < 1e-7
            end

            # Eigensolver
            D2 = gees!(A)
            D2 = eigenvalue_sorted(diagm(0 => D2))

            D1 = eigenvalue_sorted(A)

            @show norm(D0 - D1)
            @show norm(D1 - D2)

            @assert norm(D0 - D2) < 1e-7
            @assert norm(D1 - D2) < 1e-7

            # Save error for plotting
            err[itest] = norm(D0 - D2)

            itest += 1
        end
    end

    return err # For plotting
end
