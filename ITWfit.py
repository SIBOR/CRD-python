import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

file = open('D:\\CRDdata\\ITWCH4.txt')
fileOut = 'D:\\CRDdata\\wavelengthFit'

rejectBounds = [1651,1655]
fitTolerance = 0.02

T, I, W = [], [], []
for line in file:
    ln = str.split(line,'_')
    wl = float(ln[3])
    if((wl>rejectBounds[0]) & (wl<rejectBounds[1])):
        T.append(float(ln[1]))
        I.append(float(ln[5]))
        W.append(wl)

file.close()

##plane fit to wavelength as function of temperature and current
A = np.column_stack((np.ones(len(T)), I, T))
c, resid,rank,sigma = np.linalg.lstsq(A,W)

def wavelength(temp,current):
    return c[0]+c[2]*temp+c[1]*current

##Use fit to toss out outliers and refit
Tnew, Inew, Wnew = [],[],[]
for i in range(0,len(W)):
    if(abs(wavelength(T[i],I[i]) - W[i]) < fitTolerance):
        Tnew.append(T[i])
        Inew.append(I[i])
        Wnew.append(W[i])

A = np.column_stack((np.ones(len(Tnew)), Inew, Tnew))
c, resid,rank,sigma = np.linalg.lstsq(A,Wnew)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(np.array(Inew),np.array(Tnew),np.array(Wnew))
plt.show()

print(c)

##Write Coefficients to file
file = open(fileOut,'w')
file.write(str(c[0]) + '\t' + str(c[1]) + '\t' + str(c[2]))
file.close()