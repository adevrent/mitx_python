################################################################################
# CSE.0002x
# Problem Set 2: mylinsolver
# edX Username: adevrent

"""
This Python library implements Gaussian elimination to solve Ku=f, where:
    K: 2D NumPy array that is an mxm square matrix
    f: NumPy array of length m
    u: is calculated and will be a NumPy array of length m
"""

import numpy as np


def myGE(K, f):
    """
    Gaussian elimination (without vectorization)
    """
    m, n = K.shape
    assert m == n, "Non-square matrix"

    u = np.zeros(n)

    # Extended matrix with the right-hand side as last column
    A = np.zeros((n, n+1))
    A[:, 0:n] = K
    A[:, n] = f

    # Elimination of nonzero elements below the diagonal
    for i in range(n):
        assert A[i,i] != 0.0, "Zero pivot detected"
        for j in range(i+1, n):
            Lji = A[j,i] / A[i,i]
            for k in range(i+1, n+1):
                A[j,k] -= Lji*A[i,k]

    # Back substitution
    u[n-1] = A[n-1,n] / A[n-1,n-1]
    for i in range(n-2, -1, -1):
        u[i] = A[i,n]
        for j in range(i+1, n):
            u[i] -= A[i,j]*u[j]
        u[i] /= A[i,i]

    return u


def myGE_vec(K, f):
    """
    Gaussian elimination (with vectorization)
    """
    m, n = K.shape
    assert m == n, "Non-square matrix"

    u = np.zeros(n)

    # Extended matrix with the right-hand side as last column
    A = np.zeros((n, n+1))
    A[:, 0:n] = K
    A[:, n] = f

    # Elimination of nonzero elements below the diagonal
    for i in range(n):
        assert A[i,i] != 0.0, "Zero pivot detected"

        # Calculate Li which is the vector of all Lji values (for all values of j)
        #### BEGIN SOLUTION #####
        Li = A[i+1:, i] / A[i, i]  # multiplier for row i. We need to update 
                                # A[j, :] = A[j, :] - A[i, :] * Li[j] for all rows j >= i+1
        # Li is COLUMN vector!
        #### END SOLUTION ####

        # Update A[i+1:, i+1:]
        # Hint: consider using an outer product of Li and the relevant
        #       portion of the A matrix (np.outer is the NumPy function)
        #### BEGIN SOLUTION #####
        A[i+1:n, i:] -= np.outer(Li, A[i, i:])
        #### END SOLUTION ####

    # Back substitution
    u[n-1] = A[n-1,n] / A[n-1,n-1]
    for i in range(n-2, -1, -1):
        # Calculate u[i]
        # Hint: the np.dot function will be useful
        #### BEGIN SOLUTION #####
        u[i] = (A[n-2, n] - A[n-2, n-1]*u[i+1]) / A[n-2, n-2]
        #### END SOLUTION #####

    return u