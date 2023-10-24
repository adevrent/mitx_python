import numpy as np
import matplotlib.pyplot as plt

def biuniform(xmode0, xmode1, Pmode0=0.5, N=1, seed=None):
    """
    Returns an array of size N of random numbers x which have two modes 
    (i.e. the distribution is bimodal) from:
        xmode0[0] < x < xmode0[1]
        xmode1[0] < x < xmode1[1]
        
    where the probability of being in mode0 is Pmode0, and within either mode, 
    the distributions are uniform.
    
    Args:
        xmode0 (tuple): range of mode0 = xmode0[0] to xmode0[1]
        xmode1 (tuple): range of mode1 = xmode1[0] to xmode1[1]
        Pmode0 (float): probability of being in mode0
        N (shape): either an integer or shape tuple giving the size of the
                   random number array to generate
        seed (int or None): seed to be used for random number generator 
    """
    rng = np.random.default_rng(seed)
    Pmode = rng.uniform(0,1,N)
    mode0 = rng.uniform(xmode0[0],xmode0[1],N)
    mode1 = rng.uniform(xmode1[0],xmode1[1],N)
    x = np.where(Pmode>Pmode0, mode1, mode0)
 
    return x



xA = biuniform( (-2,-1), (1,2))
xB = biuniform( (-2,-1), (1,2), Pmode0=0.5)
xC = biuniform( (-2,-1), (1,2), Pmode0=0.5, N=100)
xD = biuniform( (-2,-1), (1,2), N=100)
xE = biuniform( (-2,-1), (1,2), Pmode0=0.5, N=100, seed=314)
xF = biuniform( (-2,-1), (1,2), N=100, seed=314)
