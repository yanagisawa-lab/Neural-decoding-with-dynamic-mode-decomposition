def ProjectionKernel(a, b):
    # This function calculates the projection kernel (Jihun Hamm & Daniel D. Lee, 2008).
    # The projection kernel K_p is described as,
    # k_P(Y_1, Y_2) = ||Y_1' Y_2||_F^2
    # where Y_1 and Y_2 are orthonormal matrices and Y_1' is a conjugate transpose of Y_1.
    
    A = np.conjugate(a.T)
    k_P = np.linalg.norm(A.dot(b), 'fro')
    
    return k_P