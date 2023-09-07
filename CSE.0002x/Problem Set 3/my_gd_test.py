import gd
import gd_tests
import numpy as np


""" xhist, Jhist = gd.gradient_descent(Jfunc = gd_tests.J0, xstart = np.array([0.49986708]), alpha = 0.1, nStop = 40, pdict = None)
print(xhist)
print("xhist shape =", xhist.shape) """

xhist, Jhist = gd.gradient_descent(Jfunc = gd_tests.J0, xstart = np.array([0.49993354]), alpha = 0.1, nStop = 40, pdict = None)
print(xhist)
print("xhist[0,] =", xhist[0,])

""" xhist, Jhist = gd.gradient_descent(Jfunc = gd_tests.J1, xstart = np.array([ 0.41725561, -0.20862781]), alpha = 0.1, nStop = 10, pdict = {'user': np.array([ 0.5 , -0.25])})
print(xhist)
print("xhist[0, 0] =", xhist[0, 0])
print("xhist[0, 1] =", xhist[0, 1])
print("xhist[1, 0] =", xhist[1, 0]) """

""" xhist, Jhist = gd.gradient_descent(Jfunc = gd_tests.J1, xstart = np.array([ 0.41725561, -0.20862781]), alpha = 0.1, nStop = 10, pdict = {'user': np.array([ 0.5 , -0.25])})
print(xhist)
print("xhist[0, 0] =", xhist[0, 0])
print("xhist[0, 1] =", xhist[0, 1])
print("xhist[1, 0] =", xhist[1, 0])
print("xhist[1, 1] =", xhist[1, 1])
print("xhist[2, 0] =", xhist[2, 0])
print("xhist[2, 1] =", xhist[2, 1])
print("xhist shape =", xhist.shape) """