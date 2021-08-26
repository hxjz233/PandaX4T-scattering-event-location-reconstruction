# -*- coding: utf-8 -*-

import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from scipy.stats import multivariate_normal
from skimage import restoration
from PIL import Image
from pylab import *
from sklearn.mixture import BayesianGaussianMixture

def COGcut(inputdata, cutratio):
    COGdata = inputdata.copy()
    COGx = 0
    COGy = 0
    tot = 0
    COGdata /= np.nanmax(COGdata)
    ver, hor = COGdata.shape
    for j in range(ver):
        for k in range(hor):
            if np.isnan(COGdata[j,k]) == False and COGdata[j,k] > cutratio:
                COGx += COGdata[j,k] * k
                COGy += COGdata[j,k] * j
                tot += COGdata[j,k]
    COGx /= tot
    COGy /= tot
    return COGx, COGy
    
# PandaX parameters & initialization
PMTDiameter = 3 * 25.4          #3in, to mm
PMTarea = math.pi * PMTDiameter ** 2 / 4
TPCRadius = 600

PSFsize = 200
PSFwholesize = 290

dataInProcess = []
multithreshold = 0.1

template = np.zeros([2 * PSFwholesize + 1, 2 * PSFwholesize + 1])
templatenum = 0

# PSF initialization
psf = pd.read_csv('100-1100template.csv')
psf = psf.values
PSFwhole = pd.read_csv('100-1100templatewhole.csv')
PSFwhole = PSFwhole.values
PSFtemp = pd.read_csv('resulttemp.csv')
PSFtemp = PSFtemp.values

# input data
dir = 'all'                # data directory
for file in os.listdir(path = dir):
    dataInProcess.append(file)
for i in range(len(dataInProcess)):
    data = pd.read_table(dir + '/' + dataInProcess[i], header = None, names = ['PMTID','x','y','Intensity'])
    data = data.dropna()            # Last 8 rows are not raw data, and contain NaNs
    
    x = data['x'].values.tolist()
    y = data['y'].values.tolist()
    coord = data.loc[:,['x','y']].values.tolist()
    intens = data['Intensity'].values.tolist()
    
# interpolation
    grid_x, grid_y = np.mgrid[-600:600:1, -600:600:1]
    interp = interpolate.griddata(coord, intens, (grid_x, grid_y), method='cubic')
    interp = interp.T

# preprocess
    preprocess = interp.copy()
    preprocess = np.where(np.isnan(preprocess), 1e-9, preprocess)
    preprocess = np.where(preprocess <= 0, 1e-9, preprocess)        # clear out nan/0/negatives

# estimate approximate position: COG
    intensmaxloc = np.nanargmax(preprocess)
    height, width = preprocess.shape
    maxlocy = int(intensmaxloc / height) - 600
    maxlocx = intensmaxloc % height - 600
    print("max:", maxlocx, maxlocy)
    if maxlocx ** 2 + maxlocy ** 2 < 300 ** 2:
        '''
        COGdataslice = preprocess[max(0, maxlocy + 600 - 300) : min(1200, maxlocy + 600 + 300 + 1), max(0, maxlocx + 600 - 300) : min(1200, maxlocx + 600 + 300 + 1)]
        appx, appy = COGcut(COGdataslice, 0.1)
        
        appx = int(appx)
        appy = int(appy)
        appy += max(0, maxlocy + 600 - 300)
        appy -= 600
        appx += max(0, maxlocx + 600 - 300)
        appx -= 600
        print("approximate:", appx, appy)
        '''
        finalpsf = psf.copy()
    
    # deconvolution
        preprocess /= (np.max(preprocess) * 1000)                 # compress values, prevent oversize values in deconvolution process (<= 1.)
        
        decon = restoration.richardson_lucy(preprocess, finalpsf, 40)
        '''
        print(np.where(decon<=0))
        '''
        intensmaxloc = np.nanargmax(decon)
        height, width = interp.shape
        maxlocy = int(intensmaxloc / height) - 600
        maxlocx = intensmaxloc % height - 600
        print("deconvoluted max:", maxlocx, maxlocy)
        
    # clustering after deconvolution, estimates whether is a single / multiple scattering event
        dots = np.empty([0,2])
        pointweight = decon / np.max(decon) * 20
        
        for j in range(max(0, maxlocy + 600 - 200), min(1200, maxlocy + 600 + 200 + 1)):
            for k in range(max(0, maxlocx + 600 - 200), min(1200, maxlocx + 600 + 200 + 1)):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 600, j - 600]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
        #ranheight, ranwidth = dots.shape
        #dots += np.random.uniform(-0.5, 0.5, (ranheight, ranwidth))
        bgm = BayesianGaussianMixture(n_components=2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        meansdec = bgm.means_
        weightsdec = bgm.weights_
        print(meansdec, weightsdec)
        
        # only work for 2 centers!
        if weightsdec[0] < multithreshold or weightsdec[1] < multithreshold:
            print(dataInProcess[i])
            singletemp = interp[np.arange(int(maxlocy) + 600 - PSFwholesize, int(maxlocy) + 600 + PSFwholesize + 1)][:,np.arange(int(maxlocx) + 600 - PSFwholesize,int(maxlocx) + 600 + PSFwholesize + 1)]
            #singletemp = interp[np.arange(int(yTMs) + 400, int(yTMs) + 801)][:,np.arange(int(xTMs) + 400,int(xTMs) + 801)]
            singletemp[np.isnan(singletemp)] = 0
            template += singletemp / sum(sum(singletemp))
            templatenum += 1
            
            templatesave = template.copy()
            templatesave /= templatenum
            templatesave = pd.DataFrame(templatesave)
            templatesave.to_csv("templatewhole.csv", index=False)
            print(str(templatenum) + " adopted!")
    print("finished! /" + str(i + 1))
template /= templatenum
plt.imshow(template, cmap = 'gray')
plt.show()

template = pd.DataFrame(template)
template.to_csv("templatewhole.csv", index=False)
