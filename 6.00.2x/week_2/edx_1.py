# import random
# def genEven():
#     '''
#     Returns a random even number x, where 0 <= x < 100
#     '''
#     return 2 * random.randint(0, 49)

# evenNum = genEven()
# print(evenNum)



# import random
# def deterministicNumber():
#     '''
#     Deterministically generates and returns an even number between 9 and 21
#     '''
#     return 10


import random
def stochasticNumber():
    '''
    Stochastically generates and returns a uniformly distributed even number between 9 and 21
    '''
    return 2 * random.randint(5, 10)

ans = stochasticNumber()
print(ans)