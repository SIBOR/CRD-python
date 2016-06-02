import sys, string
import numpy as np
import matplotlib
matplotlib.use('TKAgg');
import matplotlib.pyplot as plt
import random

###---------Default File name------------
#filename = '/media/win/Users/jaymz/Desktop/new_data/text.txt'

filename = None #Will prompt user if set = None

###-----------Wavelength Fit Coefficents------------
wavelengthFitFile = '/media/win/Users/jaymz/Desktop/new_data/wavelengthFit'

###-------------Parameters------------
chiSquareBounds = [1E-5, 5.9E-5]
ringdownBounds = [1.0E-5, 2.0E-4]
threeTauBounds = [1.0E-6, 1.0E-4]
wlBounds = [1598.0, 1602.0] 
maxErr = 1.0E4
correctPredictedWL = False
wavelengthTolerance = 0.1


###--------------Scan command line arguments------------------------
for i in range(0,len(sys.argv)):
    if((str.lower(sys.argv[i]) == '-i' or str.lower(sys.argv[i]) == '--input') and i < len(sys.argv) - 1):
        filename = sys.argv[i+1]
    if(str.lower(sys.argv[i]) == '-h' or str.lower(sys.argv[i]) == '--help'):
        print("\nScript for initial processing of Labview acquisition program output files.\nTakes the median of the ringdown traces in the file to produce a three collumn file of wavelength, decay constant, and error\n\nUsage:\n-i, --input\t[Input File]\t\tLocation of input file.\n--help\t\t\t\t\tDisplay this help")
        sys.exit()
        
###---------------Prompt for filename if none set----------------
def tk_file_open_dialog(prompt):
    import Tkinter as tk
    import tkFileDialog
    
    root = tk.Tk()
    root.withdraw()
    fileName = tkFileDialog.askopenfilename(title = str(prompt))
    
    if(fileName == ()):
        sys.exit()
    
    return fileName

if(filename == None):
    filename = tk_file_open_dialog("Choose Data file (output of Labview Acquisition program)")

###----------Filename for final output --------------
dataOut = filename + "_medians.csv"
dotIndex = str.rfind(filename,'.')   #peel off file extension and append _dataDump.csv to filename
if(dotIndex >=0):
    dataOut = filename[0:dotIndex] + "_medians.csv"

###Program Start

spectrum = {'wavelength'    :   [],
            'ringdownTime'  :   [],
            'error'         :   [],
            'ringdownTimeArr' : [],
            'laserInfo'     :   [],
            'allFitValues'  :   []  }


vals = []
laserInfo = []

wlCoef = []

if(correctPredictedWL):
    wlFile = open(wavelengthFitFile,'r')
    for coeff in str.split(wlFile.read(),'\t'):
        wlCoef.append(float(coeff))
    wlFile.close()
    
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
    
file = open(filename)

for line in file:
    if(line[0] != '#'):
        if(line[0:2] == '$$'):
            if(len(vals) > 0):
                val = []
                for v in vals:
                    if((v[1] > -1.0/ringdownBounds[0]) & (v[1] < -1.0/ringdownBounds[1])):
                        val.append(v[1])
                meds = bootstrapMedian(val)
                #meds = [np.mean(val), 0.0]
                if(meds[1] < maxErr):                   #reject point if bootstrap error too large
                    wlMeasure = laserInfo[0]
                    if((wlMeasure > wlBounds[0]) & (wlMeasure < wlBounds[1])):
                        spectrum['wavelength'].append( laserInfo[0] )      #wavelength measurement
                        #spectrum[0].append( wavelengthPredicted(laserInfo[1],laserInfo[3]) )      #wavelength
                        spectrum['ringdownTime'].append( meds[0] )           #median Ringdown time
                        spectrum['error'].append( meds[1] )            #error of median ringdown time
                        spectrum['ringdownTimeArr'].append( val )               #array of individual ringdown times
                        spectrum['laserInfo'].append( laserInfo )         #array of other info about acquisition
                        spectrum['allFitValues'].append( vals )             #array of all other info about individual acquisitions
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
            if((err > chiSquareBounds[0]) & (err < chiSquareBounds[1])):
                if((len(l) >= 6)):              #compatibility with datasets not containing three tau error
                    threeTau = float(l[5])
                    if((threeTau < threeTauBounds[0]) | (threeTau > threeTauBounds[1])):
                        continue
                ln = []
                for v in l:
                    ln.append(float(v))
                vals.append(ln)
file.close()

##Plane fit to wavelength based on current data
A = np.column_stack( ( np.ones(len(spectrum['wavelength'])), np.matrix(spectrum['laserInfo'])[:,2], np.matrix(spectrum['laserInfo'])[:,1] ) )
datCoef, resid,rank,sigma = np.linalg.lstsq(A,np.matrix(spectrum['laserInfo'])[:,0])
datCoef = (datCoef.flatten().tolist())[0]

##-----------Toss out wavelength measurements that deviate more than "wavelengthTolerace"---------------
if(correctPredictedWL):
    wlCoef = datCoef
    for i in range(0,len(spectrum['wavelength'])):
        wlPred = wavelengthPredicted(spectrum['laserInfo'][i][1],spectrum['laserInfo'][i][2])
        if(abs(spectrum['wavelength'][i] - wlPred) > wavelengthTolerance):
            spectrum['wavelength'][i] = wlPred


##----------Write final output file--------------   
def writeDataDumpFile(fname):    #output file format compatable with the older MATLAB scripts. This file type should be entirely obsolete now since all analysis is now done in python.
    file = open(fname,'w')
    for i in range(0,len(spectrum['wavelength'])):
        predictedWL = wavelengthPredicted(spectrum['laserInfo'][i][1], spectrum['laserInfo'][i][3])
        trig = str(spectrum['laserInfo'][i][3])
        if(len(spectrum['allFitValues'][i][0]) > 4): #compatibility with older data sets that didnt record trigger values
            trig = str(np.mean(np.matrix(spectrum['allFitValues'][i])[:,4]))
        file.write(str(spectrum['wavelength'][i]) + "," + str(-1.0/spectrum['ringdownTime'][i]) + "," + str(spectrum['laserInfo'][i][2]) + "," + "0" + "," + "0" + "," + "0" + "," + str(spectrum['laserInfo'][i][1]) + "," + str(spectrum['laserInfo'][i][3]) + "," + str(predictedWL) + "," + str(spectrum['laserInfo'][i][0]) + "," + trig + "\n")
    file.close()
    
def writeMedianFile(fname):
    file = open(fname, 'w')
    file.write("#Median ringdown times for the input data set: " + filename + "\n")
    file.write("#Note that the ringdown times are the fit coefficient, defined as -1/tau, where tau is the 1/e intesenisy decay time in sec\n")
    file.write("#Wavelength(nm)\tMedian Ringdown time(s^-1)\tError of median(s^-1)\n")
    for i in range(0,len(spectrum['wavelength'])):
        file.write(str(spectrum['wavelength'][i]) + "\t" + repr(spectrum['ringdownTime'][i]) + "\t" + repr(spectrum['error'][i]) + "\n")
    file.close()
    print("Median ringdown times written to "+ fname)

writeMedianFile(dataOut)

##--------Plot spectrum----------
#plt.plot(spectrum[0],spectrum[1])
#plt.errorbar(3.0*10**17/np.array(spectrum[0]),np.array(spectrum[1])*10**6, yerr=np.array(spectrum[2])*10**6,fmt='ro')
plt.errorbar(np.array(spectrum['wavelength']),np.array(spectrum['ringdownTime']), yerr=np.array(spectrum['error']),fmt='ro')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Ringdown Time (s^-1)')
plt.show()