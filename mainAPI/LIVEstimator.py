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
import localLinearRegression as LLR
#import matplotlib.pyplot as plot
''' 
step 1 local linear regression of p on X and pX
This step needs to be paralleled
'''
def residualize(X, ps):
    numObs = X.shape[0]
    resHat = np.zeros((numObs,1))
    hOptimal = LLR.optimalBandwidthSelection(X, ps)
    # for each observations, get prediction
    for iObs in range(numObs):
        xFit = LLR.polynomialFit(X,ps, ps[iObs], hOptimal)
        resHat[iObs] = X[iObs] - xFit[0]
    return resHat

'''
    this function is used to estimate the coefficient of the covariates, using LLR
'''
def LLRresidualEst(Y,X, ps):
    eY = residualize(Y, ps)
    numObs = X.shape[0]
    numVar = X.shape[1]
    
    
    eX = np.zeros((numObs, 2*numVar))
    for iX in range(numVar):
        eX[:, np.array([iX]) ] = residualize( X[:, np.array([iX]) ], ps)
        eX[:, np.array([iX+numVar]) ] = residualize(X[:, np.array([iX]) ]*ps, ps)
        
    beta = np.dot( np.linalg.inv(np.dot(eX.T, eX)), np.dot(eX.T, eY))

    # now get the residual of the Y
    xPs = X*np.kron(np.ones((1, numVar)), ps)
    residY = Y - np.dot(X, beta[0:2]) - np.dot(xPs, beta[2:])
    
    return residY, beta[2:]

'''
    this function is used to estimate the MTE, using LLR
'''

def LLRKpLIVOptimalBandwidth(yTilda, ps):
    hOptimal = LLR.optimalBandwidthSelection(yTilda, ps)
    return hOptimal

def LLRKpLIVMTE(yTilda, ps, hOptimal, betaDif, x, u):
    
    # test for size agree 
    if x.shape[1] != betaDif.shape[0]:
        print 'the size of X and beta does not agree'
        raise
    elif x.shape[0] != 1:
        print 'the x is not row vector'
        raise        
    elif betaDif.shape[1] != 1:
        print 'the beta is not column vector'
        raise        
    
    result = LLR.polynomialFit(yTilda,ps,u,hOptimal)
    mteHat = np.dot(x, betaDif) + result[2]
    return mteHat

'''
***********************************************************************************
 this function is used to estimate the MTE, using Polynomial
***********************************************************************************
'''


def LLRKpPolyParam(yTilda, ps, polyOrder):
    numObs = yTilda.shape[0]
    # first estimate a polynomial
    constant = np.ones((numObs,1))
    psPoly = ps
    
    '''
        yTilda = b0 + b1p + b2p^2
        iteration starts from order 2 to order k
    '''
    
    for iPoly in range(2,polyOrder+1):
        psPoly = np.concatenate((psPoly, ps**iPoly),1)
        
    regressor = np.concatenate((constant, psPoly),1)
    beta = np.dot( np.linalg.inv(np.dot(regressor.T, regressor)), np.dot(regressor.T, yTilda))
    # do not need to report constant
    return beta[1:]

    



    