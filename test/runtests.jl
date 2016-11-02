# Initialize the random number generator
rng = MersenneTwister()

# Triangular solvers

include("../trtrs.jl")

# Size of the matrix
n = 32

# Lower triangular matrix
L = zeros(Float64,n,n)
# Filling the matrix with random integer entries
for j=1:n # Column j
    L[j,j] = 1 # Should be non-zero
    L[j+1:n,j] = rand(rng, -2:2, n-j)
end

# Initializing the right-hand side
xe = rand(rng, 0:9, n) # This will be our solution
b = L * xe

x = trtrsRow(L, b)
# Let's check the result
@assert x == xe

x = trtrs(L, b)
@assert x == xe


# LU solvers
include("../getrf.jl")

# Random initialization of matrix A
L = zeros(Float64,n,n)
U = zeros(Float64,n,n)
for i=1:n
    L[i,i] = 1
    L[i+1:n,i] = rand(rng, -2:2, n-i)
    U[i,i] = rand(rng, 1:2)
    U[i,i+1:n] = rand(rng, -2:2, n-i)
end
A = L * U
A0 = copy(A)

# Initializing the right-hand side
xe = rand(rng, 0:9, n) # This will be our solution
b = A * xe

# Test our solvers
map([getrfOuter!, getrfAxpy!, getrfDot!]) do solver
    A = copy(A0)
    solver(A)
    # Solve
    x = getrs(A, b)
    @assert x == xe
end


# LU with pivoting
# Random initialization of matrix A
L = zeros(Float64,n,n)
U = zeros(Float64,n,n)
P = randperm(rng,n) # Randow row permutation
for i=1:n
    L[P[i],i] = 3 # Largest entry in the column
    L[P[i+1:n],i] = rand(rng, -2:2, n-i)
    U[i,i] = rand(rng, 1:2)
    U[i,i+1:n] = rand(rng, -2:2, n-i)
end
A = L * U
A0 = copy(A)

# Initializing the right-hand side
xe = rand(rng, 0:9, n) # This will be our solution
b = A * xe

A = copy(A0)
P = getrf!(A)
# Solve
x = getrs(A, P, b)
@assert x == xe


# Rank revealing
# Initialize matrix
d = (1.0/2.0).^(n-1:-1:0)
Q, = qr(rand(rng, n,n))
A = Q * diagm(d);

# Testing rook pivoting kernel
b = rand(rng, n)
x1 = (A1 = copy(A); P = getrf!(A1); getrs(A1, P, b))
x2 = (A1 = copy(A); (P_row, P_col) = getrfRook!(A1); getrs(A1, P_row, P_col, b))
@assert norm(x1-x2)/norm(x2) < 10*eps(Float64)


# Cholesky factorization
include("../potrf.jl")

# Random initialization of matrix A
G = zeros(Float64,n,n)
for i=1:n
    G[i,i] = rand(rng, 1:2)
    G[i+1:n,i] = rand(rng, -2:2, n-i)
end
A = G * G.'
A0 = copy(A)

# Initializing the right-hand side
xe = rand(rng, 0:9, n) # This will be our solution
b = A * xe

A = copy(A0)
potrf!(A)
# Solve
x = potrs(A, b)
@assert norm(x - xe) == 0

# Cholesky factorization
include("../geqrf.jl")

n = 32
x = rand(rng, n)
beta, v = house(x)
y = zeros(n); y[1] = norm(x)
Px = x - beta * dot(v,x) * v
if Px[1] * y[1] < 0
    y = -y
end
@assert norm(Px - y) < 10 * eps(Float64)

e1 = zeros(n); e1[1] = 1.0
x = zeros(n); x[1] = 2.0
beta, v = house(x)
@assert beta == 0.0

x[1] = -2.0
beta, v = house(x)
@assert beta == 0.0

# x[2:end] very small and x1 > 0
x = eps(Float64) * rand(rng, n)
x[1] = 1.0
beta, v = house(x)
y = zeros(n); y[1] = norm(x)
Px = x - beta * dot(v,x) * v
if Px[1] * y[1] < 0
    y = -y
end
@assert norm(Px - y) < 8.0*eps(Float64)*eps(Float64)

# x[2:end] very small and x1 < 0
x = eps(Float64) * rand(rng, n)
x[1] = -1.0
beta, v = house(x)
y = zeros(n); y[1] = norm(x)
Px = x - beta * dot(v,x) * v
if Px[1] * y[1] < 0
    y = -y
end
@assert norm(Px - y) == 0

# QR factorization

function test_QR(A)
    m = size(A,1)
    n = size(A,2)
    Q,R = qr(copy(A))
    geqrf!(A)
    for i=1:min(m,n)
        if i <= min(m,n)
            if R[i,i] * A[i,i] < 0
                R[i,:] = -R[i,:]
            end
        end
        @assert norm(R[i,i:end] - A[i,i:end]) < 1e3*eps(Float64)
    end
end

# Square matrix
n = 64
test_QR(rand(rng, n, n))

# Fat matrix
m = 32; n = 64
test_QR(rand(rng, m, n))

m = 2; n = 1024
test_QR(rand(rng, m, n))

# Thin matrix
m = 64; n = 32
test_QR(rand(rng, m, n))

m = 1024; n = 2
test_QR(rand(rng, m, n))

# Test Givens' rotations
n = 64
A = triu(rand(rng, n,n),-1)
Q,R = qr(copy(A))
for k=1:n-1
    c, s = givens(A[k,k], A[k+1,k])
    # Apply the Givens rotation to row k and k+1
    for j=k:n
        A[k,j], A[k+1,j] =
            ( c * A[k,j] - s * A[k+1,j],
              s * A[k,j] + c * A[k+1,j] )
    end
end
for i=1:n
    if R[i,i] * A[i,i] < 0
        R[i,:] = -R[i,:]
    end
    @assert norm(R[i,i:end] - A[i,i:end]) < 1e2*eps(Float64)
end


# CGS
function orthogonal_matrix(m,n)
    # Building an orthogonal matrix Q
    Q = zeros(m,n)
    for j=0:n-1, i=0:m-1
        Q[i+1,j+1] = cos(π*(2i+1)*j/2m)
    end
    for j=1:n
        Q[:,j] /= norm(Q[:,j])
    end
    return Q
end

# Testing QR factorization using CGS
n = 16
m = 32

# Building an orthogonal matrix Q
Q = orthogonal_matrix(m,n)

# Initializing an upper triangular matrix R
R = triu(Float64[ i/j for i=1:n, j=1:n ])

# Matrix A
A = Q*R
# QR factorization
RGS = geqrfCGS!(A)
# The factor Q is stored in A

# These matrices should now be equal
@assert norm(Q-A) < 1e2 * eps(Float64)
@assert norm(R-RGS) < 1e2 * eps(Float64)

# MGS

# Matrix A
A = Q*R
# QR factorization
RGS = geqrfMGS!(A)
# The factor Q is stored in A

# These matrices should now be equal
@assert norm(Q-A) < 1e2 * eps(Float64)
@assert norm(R-RGS) < 1e2 * eps(Float64)


# Test HDF5
n = 512
A = rand(rng, n, n)
Pkg.add("HDF5")
using HDF5
h5open("data.h5", "w") do file
    @write file A # alternatively, say "@write file A"
end
A0 = h5read("data.h5", "A")
@assert norm(A - A0) == 0
