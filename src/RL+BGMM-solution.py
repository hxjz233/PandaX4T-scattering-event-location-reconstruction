# -*- coding: utf-8 -*-

import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from skimage import restoration
from PIL import Image
from pylab import *
from sklearn.mixture import BayesianGaussianMixture
import time

dir = 'testdata'                # data directory
multithreshold = 0.22
PSFsize = 200
deconiter = 200

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

PSFwholesize = 290

dataInProcess = []
dataoversize = []
datafailed = []
dataasymmetric = []

grid_x, grid_y = np.mgrid[-600:600:1, -600:600:1]
grid_x_mirror, grid_y_mirror = np.mgrid[-100:1300:1, -400:400:1]

# PSF initialization
psf = pd.read_csv('templatewhole.csv')
psf = psf.values
finalpsf = psf.copy()
finalpsf = finalpsf[PSFwholesize - PSFsize: PSFwholesize + PSFsize + 1, PSFwholesize - PSFsize: PSFwholesize + PSFsize + 1]

# input data

for file in os.listdir(path = dir):
    dataInProcess.append(file)
for i in range(len(dataInProcess)):
#for i in range(2,3):
    processstarttime = time.time()
    print(dataInProcess[i])
    euler = np.array([[1,0],[0,1]])
    eventtype = ''

    data = pd.read_table(dir + '/' + dataInProcess[i], header = None, names = ['PMTID','x','y','Intensity'])
    zdata = data.iloc[-2,0]
    xCOG = eval(data.iloc[-8,0].split(" = ")[1])
    yCOG = eval(data.iloc[-7,0].split(" = ")[1])
    xTMs = eval(data.iloc[-6,0].split(" = ")[1])
    yTMs = eval(data.iloc[-5,0].split(" = ")[1])
    xPAF = eval(data.iloc[-4,0].split(" = ")[1])
    yPAF = eval(data.iloc[-3,0].split(" = ")[1])
    data = data.dropna()            # Last 8 rows are not raw data, and contain NaNs
    
    x = data['x'].values.tolist()
    y = data['y'].values.tolist()
    coord = data.loc[:,['x','y']].values.tolist()
    intens = data['Intensity'].values.tolist()
    
# interpolation    
    interp = interpolate.griddata(coord, intens, (grid_x, grid_y), method='cubic',fill_value=10)
    interp = interp.T

# preprocess
    preprocess = interp.copy()
    preprocess = np.where(np.isnan(preprocess), 1e-9, preprocess)
    preprocess = np.where(preprocess <= 0, 1e-9, preprocess)        # clear out nan/0/negatives

# estimate approximate position: COG/max
    intensmaxloc = np.nanargmax(preprocess)
    height, width = preprocess.shape
    maxlocy = int(intensmaxloc / height) - 600
    maxlocx = intensmaxloc % height - 600
    #print("max:", maxlocx, maxlocy)
    
    if maxlocx ** 2 + maxlocy ** 2 >= 560 ** 2:          # oversize event
        appx = maxlocx
        appy = maxlocy
        eventtype = 'oversized'            
        print("Event in " + dataInProcess[i] + " seemed to occur somewhere between the center of the farthest PMT and the wall. No algorithm found valid for this type of event yet. You can refer to x = " + str(appx) +", y = " + str(appy) + ", however, for a single oversized event, by 3 cm radial error\n")
        dataoversize.append(dataInProcess[i])
        #visualization skipped
        processendtime = time.time()
        print(str(i + 1) + " finished! elapsed time: " + str(processendtime - processstarttime))
        continue
    else:
        COGdataslice = preprocess[max(0, maxlocy + 600 - 300) : min(1200, maxlocy + 600 + 300 + 1), max(0, maxlocx + 600 - 300) : min(1200, maxlocx + 600 + 300 + 1)]
        appx, appy = COGcut(COGdataslice, 0.1)
        
        appx = int(appx)
        appy = int(appy)
        appy += max(0, maxlocy + 600 - 300)
        appy -= 600
        appx += max(0, maxlocx + 600 - 300)
        appx -= 600
    print("approximate:", appx, appy)
    
# estimate using margin model or central model
    if appx ** 2 + appy ** 2 < 300 ** 2:
        eventtype = 'central'
    else:
        eventtype = 'margin'
        angle = math.degrees(math.atan2(appy, appx))
            
        rotated = interp.copy()
        #mirror
        im = Image.fromarray(rotated)
        im = im.rotate(angle)
        im = im.crop((500,200,1200,1000))
        im_flipped = im.transpose(Image.FLIP_LEFT_RIGHT)
        im = array(im)
        im_flipped = array(im_flipped)
        
        rotated = np.append(im, im_flipped, axis = 1)
        coord_mirror = []
        intens_mirror = []
        for j in range(-400,400):
            for k in range(-100,1300):
                if np.isnan(rotated[j + 400,k + 100]) == False and rotated[j + 400,k + 100] > 10:
                    coord_mirror.append([j,k])
                    intens_mirror.append(rotated[j + 400,k + 100])
        #interpolation
        interp_mirror = interpolate.griddata(coord_mirror, intens_mirror, (grid_y_mirror, grid_x_mirror), method='linear')
        interp_mirror = interp_mirror.T
        
        preprocess = interp_mirror.copy()
        preprocess = np.where(np.isnan(preprocess), 1e-9, preprocess)
        preprocess = np.where(preprocess <= 0, 1e-9, preprocess)

    print(eventtype)
# deconvolution
    preprocess /= (np.max(preprocess) * 1e3)# compress values, prevent oversize values in deconvolution process (<= 1.)
    decon = restoration.richardson_lucy(preprocess, finalpsf, deconiter)
    if np.isnan(decon).any():
        print("NaNs found in deconvoluted signal. Probably caused by 0/0s in deconvolution, please try reducing iteration or examining any process that may cause abrupt change in signal value or checking if nonpositives or NaN showed up at any specific round of deconvolution")
        datafailed.append(dataInProcess[i])
        continue
# clustering after deconvolution, estimates whether is a single / multiple scattering event
    eventnum = 2
    dellist = []
    if eventtype == 'central':       
        dots = np.empty([0,2])
        pointweight = decon / np.max(decon) * 20
        
        for j in range(max(0, appy + 600 - 100), min(1200, appy + 600 + 100 + 1)):
            for k in range(max(0, appx + 600 - 100), min(1200, appx + 600 + 100 + 1)):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 600, j - 600]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
        bgm = BayesianGaussianMixture(n_components=2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        #bgm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter = 2000).fit(dots)
        meansdec = bgm.means_
        weightsdec = bgm.weights_
        #print(meansdec, weightsdec)
        
        meansdecunprocessed = np.swapaxes(meansdec, 0, 1).copy()
        weightsdecunprocessed = weightsdec.copy()
        
        for j in range(len(weightsdec)):
            if weightsdec[j] < multithreshold:
                dellist.append(j)
        meansdec = np.delete(meansdec, dellist, axis = 0)
        weightsdec = np.delete(weightsdec, dellist)
        
        eventnum = len(weightsdec)
        meansdec = np.swapaxes(meansdec, 0, 1)
        
    elif eventtype == 'margin':
        decon = decon[:,0:700]
        dots = np.empty([0,2])
        pointweight = decon / np.max(decon) * 20
        
        intensmaxloc = np.nanargmax(decon)
        maxlocy = int(intensmaxloc / 700) - 400
        maxlocx = intensmaxloc % 700 - 100
        
        #print("deconvoluted max:", maxlocx, maxlocy)
        # warning: may not evaluate events more than 10cm away from each other
        for j in range(max(0, maxlocy + 400 - 100), min(800, maxlocy + 400 + 100)):
            for k in range(max(0, maxlocx + 100 - 100), min(700, maxlocx + 100 + 100)):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 100, j - 400]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
        bgm = BayesianGaussianMixture(n_components=2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        #bgm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter = 2000).fit(dots)
        meansdec = bgm.means_
        weightsdec = bgm.weights_
        #print(meansdec, weightsdec)
        
        meansdecunprocessed = np.swapaxes(meansdec, 0, 1).copy()
        weightsdecunprocessed = weightsdec.copy()
        
        for j in range(len(weightsdec)):
            if weightsdec[j] < multithreshold:
                dellist.append(j)
        meansdec = np.delete(meansdec, dellist, axis = 0)
        weightsdec = np.delete(weightsdec, dellist)
        
        eventnum = len(weightsdec)
        meansdec = np.swapaxes(meansdec, 0, 1)
            
# clustering on intensity image, estimates event center location
    dellist = []
    if eventtype == 'central':
        dots = np.empty([0,2])
        pointweight = preprocess / np.max(preprocess) * 2
        
        for j in range(max(0, appy + 600 - 200), min(1200, appy + 600 + 200 + 1)):
            for k in range(max(0, appx + 600 - 200), min(1200, appx + 600 + 200 + 1)):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 600, j - 600]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
            
        #bgm = BayesianGaussianMixture(n_components=2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        bgm = BayesianGaussianMixture(n_components=eventnum, random_state=0, max_iter = 2000).fit(dots)
        means = bgm.means_
        weights = bgm.weights_
        #print(means, weights)
        
        for j in range(len(weights)):
            if weights[j] < multithreshold:
                dellist.append(j)
        means = np.delete(means, dellist, axis = 0)
        weights = np.delete(weights, dellist)
        means = np.swapaxes(means, 0, 1)
        
    elif eventtype == 'margin':
        dots = np.empty([0,2])
        pointweight = preprocess / np.max(preprocess) * 2
        
        for j in range(800):
            for k in range(1400):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 100, j - 400]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
            
        bgm = BayesianGaussianMixture(n_components=eventnum*2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        #bgm = BayesianGaussianMixture(n_components=eventnum*2, random_state=0, max_iter = 2000).fit(dots)
        means = bgm.means_
        weights = bgm.weights_
        #print(means, weights)
        
        for j in range(len(weights)):
            #if weights[j] < multithreshold / 2:
            #    dellist.append(j)
            if means[j, 0] ** 2 + means[j, 1] ** 2 > 600 ** 2:
                dellist.append(j)
        means = np.delete(means, dellist, axis = 0)
        weights = np.delete(weights, dellist)
        weights *= 2
        means = np.swapaxes(means, 0, 1)
        
        #rotate to original frame
        angle = math.radians(angle)
        euler = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        means = np.dot(euler, means)
    
    resultsangle = 0
    if len(weights) != eventnum:
        print("It seems that a margin event was clustered with asymmetric result.")
        dataasymmetric.append(dataInProcess[i])
    if len(weights) == 2: # output angle
        meansdec = np.swapaxes(meansdec, 0, 1)
        means = np.swapaxes(means, 0, 1)
        
        vecdecclus = meansdec[0] - meansdec[1]
        vecinstclus = means[0] - means[1]
        resultsangle = math.degrees(math.acos(np.abs(np.dot(vecdecclus, vecinstclus)) / (np.linalg.norm(vecdecclus) * np.linalg.norm(vecinstclus))))
        
        meansdec = np.swapaxes(meansdec, 0, 1)
        means = np.swapaxes(means, 0, 1)
    
# COG for each event in deconvoluted image(only one event by the time)
#    COGdataslice = decon[max(0, appy + 600 - 100) : min(1200, appy + 600 + 100 + 1), max(0, appx + 600 - 100) : #min(1200, appx + 600 + 100 + 1)]
#    RLCOGx, RLCOGy = COGcut(COGdataslice, 0.01)
#    
#    RLCOGy += max(0, appy + 600 - 100)
#    RLCOGy -= 600
#    RLCOGx += max(0, appx + 600 - 100)
#   RLCOGx -= 600
    
# visualization
    plt.gcf().set_size_inches(12, 12)
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    plt.subplot(2,2,1)
    plt.imshow(interp, extent=(-600, 600, -600, 600), origin = 'lower', cmap = 'Blues')
    plt.colorbar()
    
    plt.scatter(x, y, s = 20, alpha = 0.8, c = [intens], cmap = 'gist_heat', label = 'Raw PMT intensity')
    
    plt.legend(loc = "lower left")
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Interpolated Intensity')
    
    plt.text(x = 0, y = -800, s = zdata, ha = 'center')
    
    plt.subplot(2,2,2)
    if eventtype == 'central':
        plt.imshow(decon, extent=(-600, 600, -600, 600), origin = 'lower', cmap = 'Blues', vmin=decon.min(), vmax=decon.max())
        plt.scatter(x, y, s = 20, alpha = 0.8, c = [intens], cmap = 'gist_heat')
    elif eventtype == 'margin':
        plt.imshow(decon, extent=(-100, 600, -400, 400), origin = 'lower', cmap = 'Blues', vmin=decon.min(), vmax=decon.max())
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Deconvoluted Signal\n(rotated for margin events. x may refer to radial distance)')
    
    plt.subplot(2,2,3)
    textstr = 'RL+BGMMweight:\n' + str(weightsdecunprocessed)
    plt.text(0.05,0.8,textstr,size=10,transform=ax3.transAxes,bbox=dict(boxstyle="round",alpha=0.6,ec=(0.5,0.5,1.),fc=(0.8,0.8,1.)))
    #plt.scatter(x = RLCOGx, y = RLCOGy, s = 300, alpha = 0.6, edgecolors = 'orange', marker = '*', label = 'Event center')
    plt.scatter(meansdecunprocessed[0], meansdecunprocessed[1], s = 100 * weightsdecunprocessed, edgecolors = 'orange', label = 'Clustering center')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('RL+BGMM Result (rotated for margin events)')
    plt.legend(loc = "lower left")
    if eventtype == 'central':
        plt.imshow(decon[appy + 600 - 200 : appy + 600 + 200, appx + 600 - 200 : appx + 600 + 200], extent=(appx - 200, appx + 200, appy - 200, appy + 200), origin = 'lower', cmap = 'Blues', vmin=decon.min(), vmax=decon.max())
    elif eventtype == 'margin':
        plt.imshow(decon[300:500,450:700], extent=(350, 600, -100, 100), origin = 'lower', cmap = 'Blues', vmin=decon.min(), vmax=decon.max())
    
    # comparison with other algorithms
    plt.subplot(2,2,4)
    textstr = 'BGMM2weight:\n' + str(weights) + '\nBGMMs angle:\n' + str(resultsangle) + '°'
    plt.text(0.05,0.8,textstr,size=10,transform=ax4.transAxes,bbox=dict(boxstyle="round",alpha=0.6,ec=(0.5,0.5,1.),fc=(0.8,0.8,1.)))
    meansdec = np.dot(euler, meansdec)
    plt.scatter(x, y, s = 400, alpha = 0.8, c = [intens], cmap = 'gist_heat')
    plt.scatter(x = xCOG, y = yCOG, s = 150, alpha = 0.6, edgecolors = 'black', marker = '^', label = 'COG')
    plt.scatter(x = xTMs, y = yTMs, s = 150, alpha = 0.6, edgecolors = 'black', marker = 'o', label = 'TMs')
    plt.scatter(x = xPAF, y = yPAF, s = 150, alpha = 0.6, edgecolors = 'black', marker = 's', label = 'PAF')
    plt.scatter(meansdec[0], meansdec[1], s = 300 * weightsdec, alpha = 0.6, edgecolors = 'orange', marker='*', label = 'RL+BGMM')
    plt.scatter(means[0], means[1], s = 300 * weights, alpha = 0.6, edgecolors = 'orange', marker = 'X', label = 'RL+BGMM2')
    #plt.scatter(x = RLCOGx, y = RLCOGy, s = 150, alpha = 0.6, edgecolors = 'black', marker = '*', label = 'RL')
    plt.legend(loc = "lower left")
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Comparison with Other Algorithms')
    plt.imshow(interp[appy + 600 - 200 : appy + 600 + 200, appx + 600 - 200 : appx + 600 + 200], extent=(appx - 200, appx + 200, appy - 200, appy + 200), origin = 'lower', cmap = 'Blues' )
    
    plt.suptitle('File Name:' + dataInProcess[i] + '\nEvent Type: ' + eventtype + '    ' + str(eventnum) + ' event(s) occured')
    
    #plt.show()
    # Save plot. Requires folder "result" before running. Replaces '.txt' in output file name.
    plt.savefig('result/multivalidity/' + dataInProcess[i][:-4] + '.png')
    
    plt.clf()
    
    # Save deconvoluted image in values
    #result = pd.DataFrame(decon)
    #result.to_csv("result.csv", index=False)
    
    processendtime = time.time()
    print(str(i + 1) + " finished! elapsed time: " + str(processendtime - processstarttime))
    
print("all finished!\nOversize events:\n" + '\n'.join(dataoversize) + "\nEvents failed deconvolution:\n" + '\n'.join(datafailed) + "\nAsymmetric clusters:\n" + '\n'.join(dataasymmetric))
