import sys, string
import numpy as np
import matplotlib.pyplot as plt
import random

###---------Default File name------------
filename = '../new_data/5-30-2016.txt'
#filename = input("Enter Input file name: ") ##User file input

###-----------Wavlength Fit Coefficents------------
wavelengthFitFile = '../wavelengthFit'

###-------------Parameters------------
chiSquareBounds = [10**-5, 5.9*10**-5]
reject_chiSquare = False

ringdownBounds = [10**-5, 30*10**-4]
reject_ringdown  = False

threeTauBounds = [10**-6, 1*10**-4]
reject_threeTau  = False

wlBounds = [1650, 1655]
reject_wlBounds  = True
 
maxErr = 7*10**-4
correctPredictedWL = False
wavelengthTolerance = 0.1

###----------Final Filename output for Matlab Fitting---------------
dataOut = filename + "_dataDump.csv"
if(filename[::-1].index('.') > filename[::-1].index('/')):
    dotIndex = filename[::-1]
    dataOut = filename[0:dotIndex] + '_dataDump.csv'

###Program Start
file = open(filename,'r')

spectrum = [[],[],[],[],[],[]]
vals = []
laserInfo = []

wlFile = open(wavelengthFitFile,'r')
wlCoef = []
for coeff in str.split(wlFile.read(),'\t'):
    wlCoef.append(float(coeff))
    
def wavelengthPredicted(temperature, current):
    return wlCoef[0]+wlCoef[2]*temperature + wlCoef[1]*current

def bootstrapMedian(dat):
    N = len(dat)
    medians = []
    mediansSqr = []
    for i in range(0,N):
        synthetic = []
        for j in range(0,N):
            synthetic.append(dat[random.randint(0,N-1)])
        medians.append(np.median(synthetic))
        mediansSqr.append([medians[i]**2])
    av=np.mean(medians)
    avSqr=np.mean(mediansSqr)
    
    return [av, np.sqrt(avSqr - av**2)]

for line in file:
    if(line[0] != '#'):
        if(line[0:2] == '$$'):
            if(len(vals) > 0):
                val = []
                for v in vals:
                    if(not reject_ringdown or ((v[1] > -1.0/ringdownBounds[0]) and (v[1] < -1.0/ringdownBounds[1]))):
                        if(v[1]!=0):
                            val.append(-1.0/v[1])
                meds = bootstrapMedian(val)
                #meds = [np.mean(val), 0.0]
                if(meds[1] < maxErr):                   #reject point if bootstrap error too large
                    wlMeasure = laserInfo[0]
                    if(not reject_wlBounds or ((wlMeasure > wlBounds[0]) & (wlMeasure < wlBounds[1]))):
                        spectrum[0].append( laserInfo[0] )      #wavelength measurement
                        #spectrum[0].append( wavelengthPredicted(laserInfo[1],laserInfo[3]) )      #wavelength
                        spectrum[1].append( meds[0] )           #median Ringdown time
                        spectrum[2].append( meds[1] )           #error of median ringdown time
                        spectrum[3].append( val )               #array of individual ringdown times
                        spectrum[4].append( laserInfo )         #array of other info about acquisition
                        spectrum[5].append( vals )             #array of all other info about individual acquisitions
            vals = []
            ln = str.split(line[2:],'_')
            laserInfo = []
            for v in ln:
                if((v != 'L') & (v != 'T') & (v != 'TR') & (v != 'I') & (v != 'IS') & (v != '')):
                    if((v != 'WavemeterDisabled') & ('LAOP' not in v) & ('Error' not in v)):
                        laserInfo.append(float(v))
                    else:
                        laserInfo.append(0.0)
        else:
            l = str.split(line,'\t')
            err = float(l[3])            
            if(not reject_chiSquare or ((err > chiSquareBounds[0]) & (err < chiSquareBounds[1]))):
                if((len(l) >= 6)):              #compatibility with datasets not containing three tau error
                    threeTau = float(l[5])
                    if(reject_threeTau and ((threeTau < threeTauBounds[0]) | (threeTau > threeTauBounds[1]))):
                        continue
                ln = []
                for v in l:
                    ln.append(float(v))
                vals.append(ln)
file.close()

##Plane fit to wavelength based on current data
A = np.column_stack( ( np.ones(len(spectrum[0])), np.matrix(spectrum[4])[:,2], np.matrix(spectrum[4])[:,1] ) )
datCoef, resid,rank,sigma = np.linalg.lstsq(A,np.matrix(spectrum[4])[:,0])
datCoef = (datCoef.flatten().tolist())[0]

##-----------Toss out wavelength measurements that deviate more than "wavelengthTolerace"---------------
if(correctPredictedWL):
    wlCoef = datCoef
    for i in range(0,len(spectrum[0])):
        wlPred = wavelengthPredicted(spectrum[4][i][1],spectrum[4][i][2])
        if(abs(spectrum[0][i] - wlPred) > wavelengthTolerance):
            spectrum[0][i] = wlPred


##----------Write final output file--------------   
file = open(dataOut,'w')
for i in range(0,len(spectrum[0])):
    predictedWL = wavelengthPredicted(spectrum[4][i][1], spectrum[4][i][3])
    trig = str(spectrum[4][i][3])
    if(len(spectrum[5][i][0]) > 4): #compatibility with older data sets that didnt record trigger values
        trig = str(np.mean(np.matrix(spectrum[5][i])[:,4]))
    file.write(str(spectrum[0][i]) + "," + str(-1.0/spectrum[1][i]) + "," + str(spectrum[4][i][2]) + "," + "0" + "," + "0" + "," + "0" + "," + str(spectrum[4][i][1]) + "," + str(spectrum[4][i][3]) + "," + str(predictedWL) + "," + str(spectrum[4][i][0]) + "," + trig + "\n")
file.close()

##--------Plot spectrum----------
#plt.plot(spectrum[0],spectrum[1])
#plt.errorbar(3.0*10**17/np.array(spectrum[0]),np.array(spectrum[1])*10**6, yerr=np.array(spectrum[2])*10**6,fmt='ro')
plt.errorbar(np.array(spectrum[0]),np.array(spectrum[1])*10**6, yerr=np.array(spectrum[2])*10**6,fmt='ro')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Ringdown Time (usec)')
plt.show()