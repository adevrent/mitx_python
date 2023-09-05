import numpy as np
import matplotlib.pyplot as plt

def GaussianWells(xv):
    """
    Multiple-well Gaussian function

    Args:
        xv (1D numpy array): x=xv[0], y=xv[1]
        
    Returns:
        f(x,y) (float): value of the function at (x,y)
    """
    # location of the wells
    xW = np.array([[0.3, 0.25], [0.25, 0.8], [0.75, 0.75]])
    # steepness of the wells (smaller means steeper)
    sigma = np.array([0.2,0.2,0.2])
    # amplitudes of the wells
    A = np.array([1,1,1.2])
    f = 0.
    for i in np.arange(len(xW)):
        f -= A[i] * np.exp(-(np.linalg.norm(xv - xW[i,:]))**2 / (2*sigma[i]**2))
    return f


if __name__ == '__main__':
    xx = np.linspace(0.,1.,101)
    yy = np.linspace(0.,1.,101)
    f = np.zeros((len(yy),len(xx)))
    
    # Loop over all (x,y) points and evaluate objective function
    xv = np.zeros(2)
    for j in range(len(xx)):
        for i in range(len(yy)):
            xv[0] = xx[j]
            xv[1] = yy[i]
            f[i,j] = GaussianWells(xv)
            
    # Plot contours.  Include contour labels.  
    fig, axs = plt.subplots()
    Cf = axs.contour(xx,yy,f)
    axs.clabel(Cf)
    axs.set_xlabel('x')
    axs.set_ylabel('y')
    axs.grid(True)
    axs.axis('square')
    plt.show()

