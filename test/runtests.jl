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
Q, = qr(rand(n,n))
A = Q * diagm(d);

# Testing rook pivoting kernel
b = rand(n)
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
