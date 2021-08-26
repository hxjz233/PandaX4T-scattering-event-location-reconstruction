# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

PSFsize = 290

dataInProcess = []
grid_x, grid_y = np.mgrid[-600:600:1, -600:600:1]
template = np.zeros([2 * PSFsize + 1, 2 * PSFsize + 1])
templatenum = 0
zofinterest = [100, 1100]

dir = 'RL'                #data directory
for file in os.listdir(path = dir):
    dataInProcess.append(file)
for i in range(len(dataInProcess)):
    data = pd.read_table(dir + '/' + dataInProcess[i], header = None, names = ['PMTID','x','y','Intensity'])
    zdata = data.iloc[-2,0]
    zdata = eval(zdata.split(" = ")[1])
    xTMs = data.iloc[-6,0]
    xTMs = eval(xTMs.split(" = ")[1])
    yTMs = data.iloc[-5,0]
    yTMs = eval(yTMs.split(" = ")[1])
    if zdata > zofinterest[0] and zdata <= zofinterest[1]:
        imax = data['Intensity'].argmax()
        [[xmax,ymax]] = data.iloc[[imax],[1,2]].values
        if xmax ** 2 + ymax ** 2 < 300 ** 2:
            data = data.dropna()            #Last 8 rows are not raw data, and contain NaNs
            x = data['x'].values.tolist()
            y = data['y'].values.tolist()
            coord = data.loc[:,['x','y']].values.tolist()
            intens = data['Intensity'].values.tolist()
            
            #interpolation
            interp = interpolate.griddata(coord, intens, (grid_x, grid_y), method='cubic')
            interp = interp.T
            
            #COG
            COGx = 0
            COGy = 0
            tot = 0
            for j in range(int(ymax) + 600 - PSFsize, int(ymax) + 600 + PSFsize + 1):
                for k in range(int(xmax) + 600 - PSFsize, int(xmax) + 600 + PSFsize + 1):
                    if np.isnan(interp[j,k]) == False and interp[j,k] > 100:
                        COGx += interp[j,k] * k
                        COGy += interp[j,k] * j
                        tot += interp[j,k]
            COGx /= tot
            COGy /= tot             #array index, not coordinates
            
            #changeframe
            singletemp = interp[np.arange(int(COGy) - PSFsize, int(COGy) + PSFsize + 1)][:,np.arange(int(COGx) - PSFsize,int(COGx) + PSFsize + 1)]
            #singletemp = interp[np.arange(int(yTMs) + 400, int(yTMs) + 801)][:,np.arange(int(xTMs) + 400,int(xTMs) + 801)]
            singletemp[np.isnan(singletemp)] = 0
            template += singletemp / sum(sum(singletemp))
            templatenum += 1

template /= templatenum
plt.imshow(template, cmap = 'gray')
plt.show()

template = pd.DataFrame(template)
#template.to_csv(str(zofinterest[0]) +"-" + str(zofinterest[1]) + "template.csv", index=False)
template.to_csv(str(zofinterest[0]) +"-" + str(zofinterest[1]) + "templatewhole.csv", index=False)
