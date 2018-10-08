function getrfOuter!(A)
    n = size(A,1)
    for k=1:n
        @assert A[k,k] != 0
        for i=k+1:n
            A[i,k] /= A[k,k]
        end

        # Outer-product of column k and row k
        for i=k+1:n
            for j=k+1:n
                A[i,j] -= A[i,k] * A[k,j]
            end
        end
    end
end

function getrfAxpy!(A)
    n = size(A,1)
    for j=1:n
        for k=1:j-1
            # Linear combination of column j and k
            for i=k+1:n
                A[i,j] -= A[i,k] * A[k,j]
            end
        end

        @assert A[j,j] != 0
        for i=j+1:n
            A[i,j] /= A[j,j]
        end
    end
end

function getrfDot!(A)
    n = size(A,1)
    for j=1:n
        # Update above the diagonal
        for i=2:j
            # Dot product
            for k=1:i-1
                A[i,j] -= A[i,k] * A[k,j]
            end
        end

        # Update below the diagonal
        @assert A[j,j] != 0
        for i=j+1:n
            # Dot product
            for k=1:j-1
                A[i,j] -= A[i,k] * A[k,j]
            end
            A[i,j] /= A[j,j]
        end
    end
end

function getrf!(A)
    # Uses pivoting
    n = size(A,1)
    P = collect(1:n)
    for k=1:n
        # Find pivot
        imx = k - 1 + argmax( abs.(A[k:end,k]) ) # row with largest entry
        # Swap rows
        for j=1:n
            A[k,j],A[imx,j] = A[imx,j],A[k,j]
        end
        P[[k,imx]] = P[[imx,k]]

        # Proceed with factorization
        for i=k+1:n
            A[i,k] /= A[k,k]
        end

        for j=k+1:n, i=k+1:n
            A[i,j] -= A[i,k] * A[k,j]
        end
    end
    return P
end

function getrfRook!(A)
    # Rook pivoting for rank-revealing LU
    n = size(A,1)
    P_row = collect(1:n); P_col = collect(1:n);
    for k=1:n
        # Rook pivoting
        row = 1; row0 = 0; col = 1; col0 = 0
        while row != row0 || col != col0
            row0, col0 = row, col # Save old values
            row_A = abs.(A[row+k-1, k:end]) # Search in pivots' row
            col = argmax(row_A)
            col_A = abs.(A[k:end, col+k-1]) # Search in pivot's column
            row = argmax(col_A)
        end
        # If we reach this line, this means that the pivot is the largest
        # in its row and column.
        row += k-1; col += k-1
        # Swap rows and columns
        P_row[k], P_row[row] = P_row[row], P_row[k]
        P_col[k], P_col[col] = P_col[col], P_col[k]
        for j=1:n
            A[k,j],A[row,j] = A[row,j],A[k,j]
        end
        for i=1:n
            A[i,k],A[i,col] = A[i,col],A[i,k]
        end

        # Perform usual LU step
        if A[k,k] != 0
            for i=k+1:n
                A[i,k] /= A[k,k]
            end
            for j=k+1:n, i=k+1:n
                A[i,j] -= A[i,k] * A[k,j]
            end
        end
    end
    return P_row, P_col
end

function getrs!(A, x)
    n = length(b)
    # Solve using matrix L
    for j = 1:n
        for i = j+1:n
            x[i] -= A[i,j] * x[j]
        end
    end
    # Solve using matrix U
    for j = n:-1:1
        x[j] /= A[j,j]
        for i = 1:j-1
            x[i] -= A[i,j] * x[j]
        end
    end
end

function getrs(A, b)
    n = length(b)
    x = copy(b)
    getrs!(A, x)
    return x
end

function getrs(A, P, b)
    n = length(b)
    @assert length(P) == n
    x = Vector{Float64}(undef,n)
    for i = 1:n
        x[i] = b[P[i]]
    end
    getrs!(A, x)
    return x
end

function getrs(A, P_row, P_col, b)
    n = length(b)
    @assert length(P_row) == n && length(P_col) == n
    y = Vector{Float64}(undef,n)
    for i = 1:n
        y[i] = b[P_row[i]]
    end
    getrs!(A, y)
    x = Vector{Float64}(undef,n)
    for i = 1:n
        x[P_col[i]] = y[i]
    end
    return x
end
