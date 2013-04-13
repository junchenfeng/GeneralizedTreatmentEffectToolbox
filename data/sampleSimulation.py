''' 
---------------------------------------------------------------------------

    Copyright 2013    Junchen Feng
    
    This file is part of the Generalized Treatment Effect Toolbox. 
    
    The Generalized Treatment Effect Toolbox is free software: you can redistribute it 
    and/or modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.    
 
---------------------------------------------------------------------------

 '''


import numpy as np
import os
import csv
#import matplotlib.pyplot as plot

def sim():
    a = 0.67
    b = 0.2
    c = 1.5
    numObs = 10
    sigma = np.array([[1,-0.9],[-0.9,1]], ndmin = 2, dtype = 'float' )
    meanVec = np.array([0,0], ndmin = 1, dtype = 'float')
    
    ''' Now generate U1 and U2 '''
    U = np.random.multivariate_normal(meanVec, sigma, numObs)
    
    ''' verify the simulation | Reasonable accuracy'''
    #simMean = np.mean(U, axis = 0)
    #simCov = np.cov(U.T)
    
    #print simMean
    #print simCov
    
    ''' Now simulate Y1 and Y0 '''
    Y1 = a + b + U[:,0]
    Y0 = a + U[:,1]
    DLatent = Y1 - Y0- c
    
    D = np.array(DLatent>0, dtype = int)
    Y = Y1*D+Y0*(1-D)
    
    return Y,D

''' Now verify | Correct'''
# print Y1, Y0, DLatent, D


def underiv_data(numRowsRead):
    dataList = []
    workDir =  os.getcwd()
    f = open(workDir + '/data/underiv_data.csv' ,'rb')
    
    readed = csv.reader(f)
    countNum = 0
    for row in readed:
        countNum += 1
        dataList.append(row)
        if countNum>=numRowsRead:
            break
    
    data = np.array(dataList)
    
    
    Y = np.array(data[:,0:1], dtype = float)
    D = np.array(data[:,1:2], dtype = float)
    X = np.array(data[:,2:4], dtype = float)
    Z = np.array(data[:,4:6], dtype = float)
    
    return Y,D,X,Z
    
    

    
    