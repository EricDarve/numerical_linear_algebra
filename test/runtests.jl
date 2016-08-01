# Initialize the random number generator
rng = MersenneTwister();

@assert 1 == 1

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

# Load our triangular solvers
include("../trtrs.jl")

x = trtrsRow(L, b)
# Let's check the result
@assert x == xe

x = trtrs(L, b)
@assert x == xe
