import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pdb
import pandas as pd
import statsmodels.api as sm


#attention, strong trends make everything correlate!!


t1 = np.random.random_integers(0, high=20, size=30)
t2 = np.random.random_integers(0, high=100, size=30)


trend1 = np.linspace(-20,20.1,30)
trend2 = np.linspace(-70,70.1,30)

f = plt.figure()
ax = f.add_subplot(2,1,1)
plt.plot(t1)
ax = f.add_subplot(2,1,2)
plt.plot(t2)
plt.title('No trend')
plt.show()

f = plt.figure()
ax = f.add_subplot(2,1,1)
plt.plot(t1+trend1)
ax = f.add_subplot(2,1,2)
plt.plot(t2+trend1)
plt.title('With trend')
plt.show()


f = plt.figure()
ax = f.add_subplot(2,1,1)
plt.plot(t1+trend1)
ax = f.add_subplot(2,1,2)
plt.plot(t1+trend1)
plt.plot(t2+trend2)
plt.title('With trend')
plt.show()

print('no trend', stats.pearsonr(t1,t2))

print('Trend pearson', stats.pearsonr(t1+trend1,t2+(trend2)))
print('Trend spearman', stats.spearmanr(t1+trend1,t2+(trend2)))
print('Trend Ktau', stats.kendalltau(t1+trend1,t2+(trend2)))



f = plt.figure()
ax = f.add_subplot(1,1,1)
plt.plot(t1)
plt.plot(t2)


print('no jump', stats.pearsonr(t1,t2))
t1[0:15] = t1[0:15]-5.
t1[15::] = t1[15::]+5.

t2[0:15] = t2[0:15]-20
t2[15::] = t2[15::]+20
print('jump', stats.pearsonr(t1,t2))



