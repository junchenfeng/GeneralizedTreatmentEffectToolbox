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

import LIVEstimator as liv

import sampleSimulation as ss
import pysal.spreg as foreignModel

# get the first 100 observations for test
Y,D,X,Z = ss.underiv_data(200)

# first needs to estimate the propensity score
psResult = foreignModel.Probit(D,Z)
psHat = psResult.predy

# then get the parameter
residY, beta = liv.LLRresidualEst(Y, X, psHat)
assert (beta[0]-0.3139)<0.0001
assert (beta[1]-0.2651)<0.0001
print 'LLRresidualEst LIV passes unit test'

'''
    If the beta can be correctly estimated, then the rest of the procedure of LIV method does not need unit test
'''
polyOrder = 4
theta = liv.LLRKpPolyParam(residY, psHat, polyOrder)

assert (theta[0]-6.25)  < 0.01
assert (theta[1]+22.16) < 0.01
assert (theta[2]-30.43) < 0.01
assert (theta[3]+14.15) < 0.01

print 'LLRresidualEst Polynormial passes unit test'
