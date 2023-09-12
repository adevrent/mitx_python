import numpy as np

with open('/home/sscf/mitx_python/CSE.0002x/Unit 15/FE_Tinversion/dTdata.npy','rb') as f:
    dT = np.load(f)

# If dT[n]>0 then there is a temperature inversion on day n.
mean = dT.mean()
unbiased_std = dT.std(ddof=1)
unbiased_var = dT.var(ddof=1)
print("Size of array dT =", len(dT))
print("mean =", mean)
print("Unbiased Stdev =", unbiased_std)
print("Unbiased Variance =", unbiased_var)