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
A rudimentary LLR function
(1) Only accepts univariate case
(2) Only use tri-cube kernel
(3) Only use 1 degree polynomials
(4) search through only a limited interval, does not work well if the data has lots variation
'''


import numpy
import kernel

'''
This is the actual fit
'''
def polynomialFit(Y, X, x0, h):
    # use only polynomial 1
    numObs = X.shape[0]
    l = numpy.zeros((1,numObs))
    L = numpy.zeros((2,numObs))
    constant = numpy.ones((numObs,1))
    dX = numpy.mat(X - x0) 
    regressor = numpy.hstack([constant, dX])
    
    # make the weights 
    wVec = kernel.tricube(X, x0, h)
    
    # anti bug
    assert wVec.shape[1] == 1
    
        
    if wVec.sum() == 0:
        '''
        ERROR
        This is actually wrong, but not sure how to deal with it elegantly
        '''
        return 0, l, 0
    # standardize 
    wVec = wVec/wVec.sum()

        
    # sparsify the W
    sampleIdx = wVec.nonzero()[0]
    effectiveSampleNum = sampleIdx.shape[0]
    
    effectiveWVec = wVec[sampleIdx]
    regressor = regressor[sampleIdx,:]
    wMat = numpy.zeros((effectiveSampleNum, effectiveSampleNum))
    
    for iObs in range(effectiveSampleNum):
        wMat[iObs,iObs] = effectiveWVec[iObs]
    
    # calculate the coefficients
    XTW = regressor.H * wMat
    try:
        effectiveL = numpy.linalg.inv(  XTW * regressor )*XTW
        l[:,sampleIdx] = effectiveL[0,:]
        L[:,sampleIdx] = effectiveL
    except:
        '''
        When except happens, switch to nadaraya-waston estimator, which is just wVec
        '''
        l = wVec.T
        L = numpy.zeros((2,numObs))
        L[numpy.array([0]),:] = wVec.T
    # take the constant as forecast
    betaHat = numpy.dot(L,Y)

         

    return betaHat[0], l, betaHat[1]

'''
This is the bandwidth selection
'''
def LOORiskEstFast(Y,X,h):
    '''
        each LLR will return yHat and l(n*1) 
        fast LOO risk algo is sum((y_i-yHat_i)/(1-L(i,i)))^2
    '''
    
    '''
        The code is based on the assumption that it is a univariate
    '''
    if X.shape[1] != 1:
        print TypeError('X has to be uni-variate')
        raise
    
    # sort the X which will allow exploitation of the sparseness, if it exists
    xSortIdx = X[:,0].argsort()
    Y = Y[xSortIdx]
    X = X[xSortIdx]
    
    # find the number of observations
    numObs = X.shape[0]
    
    yHatVec = numpy.zeros((numObs, 1))
    L = numpy.zeros((numObs, numObs))
 
    for i in range(numObs):
        result = polynomialFit(Y,X,X[i],h)
        yHatVec[i] = result[0]
        L[i,:] = result[1]
    
    # now calculate the risk
    eHat = Y - yHatVec
    adjFactor = 1-numpy.diag(L)
    adjFactor.shape = (numObs,1)
    '''
     If using nadaraya-waston estimator, it is likely that it will result in 1-L_ii to be 0 and the eHat = 0
     One way to deal with this problem is to just pass it as NaN
     The other way is to exclude those points from risk estimation. For now, I am using this approach.
    '''
    sampleIdx = adjFactor != 0  
    
    if sampleIdx.sum() < 0.95*len(sampleIdx):
        # no observations, give it a very large error
        return 999999999999
      
    adjSquareErr = eHat[sampleIdx]/adjFactor[sampleIdx]
    adjSquareErrSqr = numpy.array(adjSquareErr, dtype = float)**2   
    
    rHat = adjSquareErrSqr.mean()
    
    if numpy.isnan(rHat) :
        print 'invalid rHat detected!'
        raise 
    
    return rHat
    
def LOORiskEstSlow(Y,X,h):
    
    raise TypeError('obsolete')
    
    
    # to make equivalent obs, sort
    xSortIdx = X[:,0].argsort()
    Y = Y[xSortIdx]
    X = X[xSortIdx]    
    
    # find the number of observations
    numObs = X.shape[0]
    
    eHat = numpy.zeros((numObs,1))
    
    # store the Yhat and l
    for i in range(numObs):
        
        # creat a temp X and Y
        testX = X
        testY = Y
        
        testX = numpy.delete(testX, i, 0)
        testY = numpy.delete(testY, i, 0)
        
        
        result = polynomialFit(testY,testX,X[i],h)
        eHat[i] = Y[i] - result[0] 
        
    rHat = numpy.mean( eHat**2 )
    return rHat



def optimalBandwidthSelection(Y, X, *args):
    '''
        for a more reasonable search grid that is robust to the outliers
        use log than transform back to level 
    '''
    # add minor value to enable 
    if len(args) == 0:
        logX = numpy.log(X-X.min()+0.00001)
        hVec = numpy.array(numpy.linspace(logX.min(), logX.max(),52))
        # cut off the first couple, which will be too small anyway
        hVec = numpy.delete(hVec, [0,1,2,3,51], 0) 
        hVec = numpy.exp(hVec)  
    elif len(args) == 1:
        hVec = args[0]
    else:
        raise TypeError('Error: hVec format is wrong.')
    
    
    numH = hVec.shape[0]
    rHatVecPP = numpy.zeros((numH,1))
   
    
    import pp, sys

    # tuple of all parallel python servers to connect with
    ppservers = ()
    
    if len(sys.argv) > 1:
        ncpus = int(sys.argv[1])
        # Creates jobserver with ncpus workers
        job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers=ppservers)

    
    # see what comes out
    jobs = [(i, job_server.submit(LOORiskEstFast, (Y, X, hVec[i]), (polynomialFit,), ("numpy","kernel",))) for i in range(numH) ]
    
    
    for i, job in jobs:
        result = job()
        rHatVecPP[i] = result

    job_server.destroy()
       
    
    #start_time = time.time()
    #for iH in range(numH):
    #    rHatVecNP[iH] = LOORiskEstFast(Y,X,hVec[iH])
    #print "NP: Time elapsed ", time.time()-start_time
    
        
    # now get the 
    optimalIdx = numpy.argmin(rHatVecPP)
    return hVec[optimalIdx]
    