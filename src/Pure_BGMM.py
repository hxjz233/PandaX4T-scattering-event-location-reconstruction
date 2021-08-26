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
import matplotlib as mpl
from scipy import linalg
import itertools

dir = 'Testoversize'                # data directory
multithreshold = 0.25
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


def plot_results(X, Y_, means, covariances, index, title):
    #splot = plt.subplot(2, 1, 1 + index)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        #v, w = linalg.eigh(covar)
        #print(v, w)
        v = [covar, covar]
        w = [[1,0],[0,1]]
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color,alpha=0.1)


    #plt.xlim(-9., 5.)
    #plt.ylim(-3., 6.)
    #plt.xticks(())
    #plt.yticks(())
    plt.title(title)


# PandaX parameters & initialization
dataInProcess = []
dataoversize = []
datafailed = []
dataasymmetric = []

grid_x, grid_y = np.mgrid[-600:600:1, -600:600:1]
grid_x_mirror, grid_y_mirror = np.mgrid[200:1000:1, -400:400:1]

color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])

# input data
for file in os.listdir(path = dir):
    dataInProcess.append(file)
for i in range(len(dataInProcess)):
#for i in range(1,2):
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
    interp = interpolate.griddata(coord, intens, (grid_x, grid_y), method='cubic')
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
    if maxlocx ** 2 + maxlocy ** 2 < 550 ** 2:
        COGdataslice = preprocess[max(0, maxlocy + 600 - 300) : min(1200, maxlocy + 600 + 300 + 1), max(0, maxlocx + 600 - 300) : min(1200, maxlocx + 600 + 300 + 1)]
        appx, appy = COGcut(COGdataslice, 0.1)
        
        appx = int(appx)
        appy = int(appy)
        appy += max(0, maxlocy + 600 - 300)
        appy -= 600
        appx += max(0, maxlocx + 600 - 300)
        appx -= 600
        print("approximate:", appx, appy)
        if appx ** 2 + appy ** 2 < 450 ** 2:
            eventtype = 'central'
        else:
            eventtype = 'margin'
    else:
        appx = maxlocx
        appy = maxlocy
        eventtype = 'oversized'
        dataoversize.append(dataInProcess[i])
        '''
        print("Event in " + dataInProcess[i] + " seemed to occur somewhere between the center of the farthest PMT and the wall. No algorithm found valid for this type of event yet. You can refer to x = " + str(appx) +", y = " + str(appy) + ", however, for a single oversized event, by 3 cm radial error\n")
        
        #visualization skipped
        processendtime = time.time()
        print(str(i + 1) + " finished! elapsed time: " + str(processendtime - processstarttime))
        continue
    '''
    print(eventtype)
    
    if eventtype != 'central':
        angle = math.degrees(math.atan2(appy, appx))
            
        rotated = interp.copy()
        #mirror
        im = Image.fromarray(rotated)
        im = im.rotate(angle)
        im = im.crop((800,200,1200,1000))
        im_flipped = im.transpose(Image.FLIP_LEFT_RIGHT)
        im = array(im)
        im_flipped = array(im_flipped)
        
        rotated = np.append(im, im_flipped, axis = 1)
        coord_mirror = []
        intens_mirror = []
        for j in range(-400,400):
            for k in range(200,1000):
                if np.isnan(rotated[j + 400,k - 200]) == False and rotated[j + 400,k - 200] > 10:
                    coord_mirror.append([j,k])
                    intens_mirror.append(rotated[j + 400,k - 200])
        #interpolation
        interp_mirror = interpolate.griddata(coord_mirror, intens_mirror, (grid_y_mirror, grid_x_mirror), method='linear')
        interp_mirror = interp_mirror.T
        
        preprocess = interp_mirror.copy()
        preprocess = np.where(np.isnan(preprocess), 1e-9, preprocess)
        preprocess = np.where(preprocess <= 0, 1e-9, preprocess)
        
# clustering on intensity image, estimates event center location
    dots = np.empty([0,2])
    pointweight = preprocess / np.max(preprocess) * 2    
    
    if eventtype == 'central':
        for j in range(max(0, appy + 600 - 200), min(1200, appy + 600 + 200 + 1)):
            for k in range(max(0, appx + 600 - 200), min(1200, appx + 600 + 200 + 1)):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k - 600, j - 600]), (int(pointweight[j,k]), 1)), axis = 0)
        #print(dots[-1:])
    
        eventnum = 2
        weights = 0
        while eventnum > 0 and np.any(weights < multithreshold):
            dpgmm = BayesianGaussianMixture(n_components=eventnum, covariance_type='spherical', max_iter = 10000,weight_concentration_prior_type='dirichlet_distribution', weight_concentration_prior=1e-3).fit(dots)
            weights = dpgmm.weights_
            print(weights)
            eventnum -= 1
    
        eventnum += 1
        means = np.swapaxes(dpgmm.means_, 0, 1)
        
    else:
        for j in range(800):
            for k in range(800):
                if pointweight[j,k] >= 1:
                    dots = np.append(dots, np.tile(np.array([k + 200, j - 400]), (int(pointweight[j,k]), 1)), axis = 0)
            #print(dots[-1:])
        
        eventnum = 2
        weights = 0
        while eventnum > 0 and np.any(weights < multithreshold / 2):
            dpgmm = BayesianGaussianMixture(n_components=eventnum*2, covariance_type='spherical', n_init=10, max_iter = 10000,weight_concentration_prior_type='dirichlet_distribution', weight_concentration_prior=1e-3).fit(dots)
            weights = dpgmm.weights_
            print(weights)
            eventnum -= 1
    
        eventnum += 1
        means = dpgmm.means_
        
        '''    
        bgm = BayesianGaussianMixture(n_components=eventnum*2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
        #bgm = BayesianGaussianMixture(n_components=eventnum*2, random_state=0, max_iter = 2000).fit(dots)
        means = bgm.means_
        weights = bgm.weights_
        #print(means, weights)
        '''
        dellist = []
        for j in range(len(weights)):
            if means[j, 0] ** 2 + means[j, 1] ** 2 > 600 ** 2:
                dellist.append(j)
        means = np.delete(means, dellist, axis = 0)
        weights = np.delete(weights, dellist)
        weights *= 2
        means = np.swapaxes(means, 0, 1)
        meansinmirror = means.copy()
        
        #rotate to original frame
        angle = math.radians(angle)
        euler = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        means = np.dot(euler, means)
    
    resultsangle = 0
    if len(weights) != eventnum:
        print("===It seems that a margin event was clustered with asymmetric result.===")
        
    #bgm = BayesianGaussianMixture(n_components=2, random_state=0, covariance_type='spherical', max_iter = 2000).fit(dots)
    
    
# visualization
    #plt.gcf().set_size_inches(6, 12)
    fig = plt.figure(figsize=(16,5.5),dpi=300)
    ax2 = fig.add_subplot(1,3,2)
    
    plt.subplot(1,3,1)
    plt.imshow(interp, extent=(-600, 600, -600, 600), origin = 'lower', cmap = 'Blues')
    plt.colorbar()
    
    plt.scatter(x, y, s = 20, alpha = 0.8, c = [intens], cmap = 'gist_heat', label = 'Raw PMT intensity')
    
    plt.legend(loc = "lower left")
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Interpolated Intensity')
    
    plt.text(x = 0, y = -800, s = zdata, ha = 'center')
    
    plt.subplot(1,3,2)
    textstr = 'BGMMweight:\n' + str(dpgmm.weights_)
    plt.text(0.05,0.9,textstr,size=10,transform=ax2.transAxes,bbox=dict(boxstyle="round",alpha=0.6,ec=(0.5,0.5,1.),fc=(0.8,0.8,1.)))
    
    if eventtype == 'central':
        plot_results(dots, dpgmm.predict(dots), dpgmm.means_, dpgmm.covariances_, 1,'Bayesian Gaussian Mixture with a Dirichlet process prior')
        plt.imshow(interp[appy + 600 - 200 : appy + 600 + 200, appx + 600 - 200 : appx + 600 + 200], extent=(appx - 200, appx + 200, appy - 200, appy + 200), origin = 'lower', cmap = 'Blues' )
    else:
        plot_results(dots, dpgmm.predict(dots), dpgmm.means_, dpgmm.covariances_, 1,'Bayesian Gaussian Mixture with a Dirichlet process prior')
        plt.imshow(interp_mirror, extent=(200, 1000, -400, 400), origin = 'lower', cmap = 'Blues' )
    
    # comparison with other algorithms
    ud = [appy + 600 - 200 , appy + 600 + 200]
    lr = [appx + 600 - 200 , appx + 600 + 200]
    if ud[0] < 0:
        ud = [0, 400]
    if ud[1] > 1200:
        ud = [800, 1200]
    if lr[0] < 0:
        lr = [0, 400]
    if lr[1] > 1200:
        lr = [800, 1200]
        
    plt.subplot(1,3,3)
    plt.scatter(x, y, s = 400, alpha = 0.8, c = [intens], cmap = 'gist_heat')
    plt.scatter(x = xCOG, y = yCOG, s = 150, alpha = 0.6, edgecolors = 'black', marker = '^', label = 'COG')
    plt.scatter(x = xTMs, y = yTMs, s = 150, alpha = 0.6, edgecolors = 'black', marker = 'o', label = 'TMs')
    plt.scatter(x = xPAF, y = yPAF, s = 150, alpha = 0.6, edgecolors = 'black', marker = 's', label = 'PAF')
    plt.scatter(means[0], means[1], s = 300 * weights, alpha = 0.6, edgecolors = 'orange', marker = 'X', label = 'BGMM')
    #plt.scatter(x = RLCOGx, y = RLCOGy, s = 150, alpha = 0.6, edgecolors = 'black', marker = '*', label = 'RL')
    plt.legend(loc = "lower left")
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Comparison with Other Algorithms')
    plt.imshow(interp[ud[0] : ud[1], lr[0] : lr[1]], extent=(lr[0] - 600, lr[1] - 600, ud[0] - 600, ud[1] - 600), origin = 'lower', cmap = 'Blues')
    
    plt.suptitle('File Name:' + dataInProcess[i] + '\nEvent Type: ' + eventtype + '    ' + str(eventnum) + ' event(s) occured')
    
    plt.show()
    # Save plot. Requires folder "result" before running. Replaces '.txt' in output file name.
    #plt.savefig('result/multivalidity/' + dataInProcess[i][:-4] + '.png')
    
    plt.clf()
    
    # Save deconvoluted image in values
    #result = pd.DataFrame(decon)
    #result.to_csv("result.csv", index=False)
    
    processendtime = time.time()
    print(str(i + 1) + " finished! elapsed time: " + str(processendtime - processstarttime))
    
print("all finished!\nOversize events:\n" + '\n'.join(dataoversize) + "\nEvents failed deconvolution:\n" + '\n'.join(datafailed) + "\nAsymmetric clusters:\n" + '\n'.join(dataasymmetric))
