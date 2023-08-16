def myGE_vec(K, f):
    """
    Gaussian elimination (with vectorization)
    """
    m,n = K.shape
    assert m == n, "Non-square matrix"

    u = np.zeros(n)

    # Extended matrix with the right-hand side as last column
    A = np.zeros((n,n+1))
    A[:,0:n] = K
    A[:,n] = f

    # Elimination of nonzero elements below the diagonal
    print("A at the beginning:")
    print(A)
    print("**************************************")
    for i in range(n):  # for each ROW i:
        assert A[i,i] != 0.0, "Zero pivot detected"

        # Calculate Li which is the vector of all Lji values (for all values of j)
        Li = A[i+1:, i] / A[i, i]  # multiplier for row i. We need to update 
                                # A[j, :] = A[j, :] - A[i, :] * Li[j] for all rows j >= i+1
        # Li is COLUMN vector!
        print("Shape of Li:", Li.shape)
        
        
        # Update A[i+1:n,i+1:]
        # Hint: consider using an outer product of Li and the relevant
        #       portion of the A matrix (np.outer is the NumPy command)
        # TODO
        A[i+1:n, i:] -= np.outer(Li, A[i, i:])
        print("A at iteration", i, ":")
        print(A.round(4))
        print("---------------------------------------------------------------")
    print("Extended matrix A' after eliminations (upper triangular):")
    print(A.round(4))

                
    # Back substitution
    u[n-1] = A[n-1,n]/A[n-1,n-1]

    for i in range(n-2,-1,-1):
        # Calculate u[i].
        # Hint: the np.dot command will be useful
        u[i] = (A[n-2, n] - A[n-2, n-1]*u[i+1]) / A[n-2, n-2]

    return u