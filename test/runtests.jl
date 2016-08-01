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
include("getrf.jl")

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
    println(string(solver) * ": PASSED")
end
