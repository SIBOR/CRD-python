from scipy import stats
pw=[]
mw=[]
for i in range(0,len(spectrum[0])):
    wm = spectrum[4][i][0]
    if((wm > 1599) & (wm < 1602)):
        mw.append(spectrum[4][i][0])
        pw.append(wavelengthPredicted(spectrum[4][i][1],spectrum[4][i][3]))
wlSlope, wlIntercept, r_value, p_value, std_err = stats.linregress(pw,mw)
plt.plot(pw,mw)
plt.plot(pw,wlSlope*np.array(pw) + wlIntercept)
plt.xlabel('predicted wavelength')
plt.ylabel('measured wavelength')
plt.show()