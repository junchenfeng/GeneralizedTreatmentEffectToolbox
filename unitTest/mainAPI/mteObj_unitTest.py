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


'''
    Mainly test for the heckit case and the the linear but non-normal case
'''

import sampleSimulation as ss
import numpy as np


Y,D,X,Z = ss.underiv_data(5000)

import mteObj as mo
x0 = np.mean(X, 0)
x0.shape = (1,2)


# calculate the MTE
testObj = mo.hackit()
testObj.assignData(Y, D, X, Z)

''' 
    test heck it
'''
linear = 1
normal = 1
LLRMTE = 0
polyOrder = 4
testObj.assignParam(linear, normal, LLRMTE, polyOrder)
testObj.estimateParam()

assert abs(testObj.betaTreat[0] - 0.24)<0.01
assert abs(testObj.betaTreat[1] - 0.8)<0.01
assert abs(testObj.betaTreat[2] - 0.4)<0.01
assert abs(testObj.betaTreat[3] + 0.012)<0.001

assert abs(testObj.betaControl[0] - 0.02)<0.01
assert abs(testObj.betaControl[1] - 0.5)<0.01
assert abs(testObj.betaControl[2] - 0.1)<0.01
assert abs(testObj.betaControl[3] - 0.05)<0.01

print 'Heck it passed unit test'

normal = 0
testObj.assignParam(linear, normal, LLRMTE, polyOrder)
testObj.estimateParam()

assert abs(testObj.bDif[0] - 0.3)<0.01
assert abs(testObj.bDif[1] - 0.3)<0.01

assert abs(testObj.theta[0] - 0.641)<0.001
assert abs(testObj.theta[1] + 1.536)<0.001
assert abs(testObj.theta[2] - 2.209)<0.001
assert abs(testObj.theta[3] + 1.119)<0.001

print 'linear polynomial passed unit test'

'''
    Notice:
    (1) MTE generating function is not tested
    (2) the nonlinear approach is unit tested in LIVEstimator
'''
