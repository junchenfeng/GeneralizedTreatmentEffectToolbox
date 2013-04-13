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

import localLinearRegression as LLR
import numpy as np

dataList = []
f = open('D:/gitRepo/grmEstimatorToolbox/mteApproach/data/testData.csv','rb')
import csv
readed = csv.reader(f)
for row in readed:
    dataList.append(row)

'''
    array - the list
'''
data = np.array(dataList, dtype = float)
Y = data[:,0:1]**(1./3.)
X = data[:,1:2]   

'''
    for unit test purpose, specify the hVec
'''
hVec = np.arange(20,100,1)
hOptimal = LLR.optimalBandwidthSelection(Y, X, hVec)
assert abs(hOptimal - 93)<1

hOptimal = LLR.optimalBandwidthSelection(Y, X)
assert abs(hOptimal - 84)<1

print('localLinearRegression passes unit test')
    
    
