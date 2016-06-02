import sys, string
import numpy as np
import matplotlib.pyplot as plt

def histogram(vals, intervals):
    maxVal, minVal = max(vals), min(vals)
    counts = [0]*(intervals + 1)
    for val in vals:
        if(abs(1.0/val) > 100):
        #print(str((float(val - minVal)*intervals/(maxVal-minVal))))
            counts[int(float(val - minVal)*intervals/(maxVal-minVal))] += 1
    xvals = []
    for i in range(0,intervals+1):
        xvals.append(minVal + float(maxVal - minVal)/intervals*(i+0.5))
    return [xvals,counts]
    
def smoothHistogram(values, intervalsMin, intervalsMax):
    hgram = [[],[]]
    vals = values
    if(type(vals).__module__ == np.__name__):
        vals = values.tolist()
    for i in range(intervalsMin, intervalsMax):
        h = histogram(vals,i)
        hgram[0]+=h[0]
        hgram[1]+=h[1]
    swapped = True
    i, n = 0, len(hgram[0])
    while swapped:
        swapped = False
        for j in range(0,n-i - 1):
            if(hgram[0][j] > hgram[0][j+1]):
                temp = [ hgram[0][j], hgram[1][j] ]
                hgram[0][j] = hgram[0][j+1]
                hgram[1][j] = hgram[1][j+1]
                hgram[0][j+1] = temp[0]
                hgram[1][j+1] = temp[1]
                swapped = True
    return hgram
    
        
for i in range(0,10):
    #hgram = smoothHistogram(spectrum[3][i],10,15)
    #plt.plot(np.array(hgram[0])*10.0**6,hgram[1])
    #counts = hgram[1]
    counts, bins, patches = plt.hist(np.array(spectrum[3][i])*10.0**6, bins = 20)
    plt.xlabel('Ringdown Time (uSec)')
    plt.ylabel('Counts')
    plt.title('Ringdown time histogram')
    med = bootstrapMedian(np.array(spectrum[3][i])*10**6)
    errHandle = plt.errorbar(med[0],max(counts)/2,xerr=med[1],fmt='ro',label = 'Median Value')
    errLegend = plt.legend(handles=[errHandle], loc=1)
    #ax = plt.gca().add_artist(errLegend)
    plt.show()