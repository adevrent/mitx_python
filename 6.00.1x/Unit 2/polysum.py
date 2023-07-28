"""
Write a function called polysum that takes 2 arguments, n and s. 
This function should sum the area and square of the perimeter of the regular polygon. 
The function returns the sum, rounded to 4 decimal places.
"""
import numpy as np

def polysum(n, s):
    """
    n : int, n > 2
    s : int or float, s > 0
    """
    area = (0.25*n*s**2) / np.tan(np.pi / n)
    perimeter = n*s
    ans = area + perimeter**2
    ans = round(ans, 4)
    return ans

def testPolysum():
    testlist_s = [1, 3.15, 28, 0.00175, 1.27895, 2.0]
    for n in range(5, 10):
        for s in testlist_s:
            print("Polysum is", polysum(n, s), "for n =", n, "and s =", s)