mutable struct SparseMatrixCSR
    m::Int # Number of rows
    n::Int # Number of columns
    rowptr::Vector{Int} # Starting index for each row
    colval::Vector{Int} # Column indices
    nzval::Vector{Float64} # Matrix entries
end

# Convert a CSC matrix to CSR format
# CSC = compressed sparse column
# CSR = compressed sparse row
function SparseMatrixCSR(A_CSC::SparseMatrixCSC)
    m = A_CSC.m
    n = A_CSC.n

    # Counting the number of entries in each row
    rowptr = zeros(Int,m+1)
    for j=1:A_CSC.n
        for k=A_CSC.colptr[j]:A_CSC.colptr[j+1]-1
            rowptr[A_CSC.rowval[k]] += 1
        end
    end

    # Beginning index for each row
    rowptr = cumsum(rowptr) .+ 1
    rowptr = vcat([1],rowptr[1:end-1])

    colval = Vector{Int}(undef,rowptr[end])
    nzval = Vector{Float64}(undef,rowptr[end])
    for j=1:A_CSC.n
        for k=A_CSC.colptr[j]:A_CSC.colptr[j+1]-1
            i = A_CSC.rowval[k]
            k_CSR = rowptr[i]
            colval[k_CSR] = j
            nzval[k_CSR]  = A_CSC.nzval[k]
            rowptr[i] += 1
        end
    end
    # Reset to correct value
    rowptr = vcat([1],rowptr[1:end-1])

    return SparseMatrixCSR(m,n,rowptr,colval,nzval)
end

"""Multiplication operator for a sparse matrix A and a dense vector x"""
function Base.:*(A::SparseMatrixCSR, x::Array{Float64,1})
    # y = A * x
    y = Vector{Float64}(undef,A.m)
    for i=1:A.m
        y[i] = 0.0
        for k=A.rowptr[i]:A.rowptr[i+1]-1
            y[i] += A.nzval[k] * x[A.colval[k]]
        end
    end
    return y
end
