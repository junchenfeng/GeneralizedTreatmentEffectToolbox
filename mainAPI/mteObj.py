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
from scipy.stats import norm
import pysal.spreg as foreignModel
import LIVEstimator as liv
    

class hackit():
    '''
    tribute to the heckit
    '''
    def __init__(self):
        self.Y = []
        self.D = []
        self.X = []
        self.Z = []
        
        self.ps = []
        self.gamma = []
        

        
        self.linear = []
        self.normal = []
        
        self.LLRMTE = []
        
        self.polyOrder = []
 
        self.betaTreat = []
        self.betaControl = [] 
 
        
        self.bDif = []
        self.theta = []
        
        self.residY = []
        self.hOptimal = []
    
    def assignData(self,Y,D,X,Z):
        
        
        # anti-bug
        # all the inputs should be numpy array
        assert(isinstance(Y, np.ndarray))
        assert(isinstance(D, np.ndarray))
        assert(isinstance(X, np.ndarray))
        assert(isinstance(Z, np.ndarray))
        # they should have same observation length
        numObs_Y = Y.shape[0]
        numObs_D = D.shape[0]
        numObs_X = X.shape[0]
        numObs_Z = Z.shape[0]
        
        if numObs_Y == numObs_D and numObs_Y == numObs_X and numObs_Y == numObs_Z:
            
        
            self.Y = Y
            self.D = D
            self.X = X
            self.Z = Z
            
            
        else:
            print 'Sizes of the data are not consistent'
            raise 
    
    
    def assignParam(self, linear, normal, LLRMTE, polyOrder):
        '''
            The structure is given in the design.pdf
        '''
        self.linear = linear
        self.normal = normal
        self.LLRMTE = LLRMTE
        self.polyOrder = polyOrder
        
        if self.linear == 1 and self.normal == 1:
            print 'heckit model'
        elif self.linear == 1 and self.normal == 0:
            print 'linear polynomial model'    
        elif self.linear == 0 and self.LLRMTE == 0:
            print 'Nonlinear polynomial model'  
        elif self.linear == 0 and self.LLRMTE == 1: 
            print 'Nonlinear LIV model' 
        else:
            raise TypeError('unknown mode!')
                 
    def estimateParam(self):
        # estimate the propensity score
        self.ps, self.gamma = self._estPscore()
        print 'propensity score estimated.'
        
        if self.linear == 1 and self.normal == 1:
            # now estimate the parameter for the outcome function
            treatedIdx = self.D.nonzero()[0]
            controlIdx = np.where(self.D == 0)[0]
            self.betaTreat = self._normal_outcomeParamEst(self.Y[treatedIdx], self.X[treatedIdx,:], self.D[treatedIdx], self.ps[treatedIdx])
            print 'Y1 estimated.'
            self.betaControl = self._normal_outcomeParamEst( self.Y[controlIdx], self.X[controlIdx,:], self.D[controlIdx], self.ps[controlIdx])
            print 'Y0 estimated.'
        elif self.linear == 1 and self.normal == 0:
            # use polynomials to estimate 
            b0, b1, self.bDif, self.theta = self._polynomialEst()
            print 'Y estimated.'
            
        elif self.linear == 0 and self.LLRMTE == 0:
            self.residY, self.bDif = liv.LLRresidualEst(self.Y, self.X, self.ps)
            print 'beta estimated.'
            self.theta = liv.LLRKpPolyParam(self.residY, self.ps, self.polyOrder)
            print 'theta estimated.'
            
        elif self.linear == 0 and self.LLRMTE == 1: 
            self.residY, self.bDif = liv.LLRresidualEst(self.Y, self.X, self.ps)
            print 'beta estimated.'
            # get the optimal bandwidth for the residY, ps LLR
            self.hOptimal = liv.LLRKpLIVOptimalBandwidth(self.residY, self.ps)   
            print 'optimal bandwidth estimated.'        
        
    def inverseMillerRatio(self,ps,D):
        invPs = norm.ppf(ps)
        
        imrD1 = -norm.pdf(invPs)/ps
        imrD0 = norm.pdf(invPs)/(1-ps)
        
        imr = imrD1*D + imrD0*(1-D)
        
        return imr
    
    def _estPscore(self):
        mResult = foreignModel.Probit(self.D,self.Z)
        return mResult.predy, mResult.betas
        
    def _normal_outcomeParamEst(self, Y, X, D, ps):
        # get the inverse miller ratio
        numObs = X.shape[0]
        imrVec = self.inverseMillerRatio(ps, D)
        constant = np.ones((numObs, 1))
        regressor = np.concatenate((constant,X, imrVec), 1)
        
        beta = np.dot( np.linalg.inv(np.dot(regressor.T, regressor)), np.dot(regressor.T, Y))
        return beta
    
    def _polynomialEst(self):
        numObs = self.X.shape[0]
        numVar = self.X.shape[1]
        
        constant = np.ones((numObs,1))
        
        xPs = self.X*np.kron(np.ones((1, numVar)), self.ps)
        
        psPoly = self.ps
        for iPoly in range(2,self.polyOrder+1):
            psPoly = np.concatenate((psPoly, self.ps**iPoly),1)
        
        regressor = np.concatenate((constant, self.X, xPs, psPoly),1)
        beta = np.dot( np.linalg.inv(np.dot(regressor.T, regressor)), np.dot(regressor.T, self.Y))
        
        # parse the beta
        b0 = beta[0]
        b = beta[1:numVar+1]
        bDif = beta[numVar+1:2*numVar+1]
        sigma = beta[2*numVar+1:]
        
        return b0, b, bDif, sigma
        
    def mtePredict(self, x, u):
        if self.linear == 1 and self.normal == 1:
            dif = self.betaTreat - self.betaControl
            consDif = dif[0]
            xDif = dif[1:-1]
            sigmaDif = dif[-1]
            
            mteHat = consDif + np.dot(x, xDif) + sigmaDif* norm.ppf(u)
        elif (self.linear == 1 and self.normal == 0 )  or (self.linear == 0 and self.LLRMTE == 0):
            pPoly = []
            '''
            Intentionally put polyOrder here. Because it needs to take 1 derivative
            '''
            for iP in range(1,self.polyOrder):
                pPoly.append((iP+1)*u**iP)
            pPoly = np.array(pPoly)
            
            mteHat = self.theta[0] + np.dot(x, self.bDif) + np.dot(pPoly, self.theta[1:])
        else:
            
            mteHat = liv.LLRKpLIVMTE(self.residY, self.ps, self.hOptimal, self.bDif, x, u)                       
        return mteHat
    
    def calcfPs(self,x, U):
        '''
         second input is not useful in this case.
         just for providing standard interface
        
        
        # assert x has to be a row vector 
        #assert x.shape[0] == 1
        #fPsHat = norm.pdf(self.gamma[0]+ np.dot(x, self.gamma[1:]) )
        '''
        
        '''
        This method is got from page 12 of Estimation of Treatment effects under Essential Heterogeneity
        access through http://jenni.uchicago.edu/underiv/documentation_2006_03_20.pdf
        '''
        fPsHat = np.mean(self.ps>U)
        return fPsHat
    
    def calcPs(self,x):
        '''
        # assert x has to be a row vector 
        assert x.shape[0] == 1
        
        pScoreHat = norm.cdf(self.gamma[0]+ np.dot(x, self.gamma[1:]) )
        '''
        
        '''
         This method is got from page 12 of Estimation of Treatment effects under Essential Heterogeneity
        access through http://jenni.uchicago.edu/underiv/documentation_2006_03_20.pdf       
        '''
        
        '''
        There is a theorem states the integral of CDF is the expected value
        '''
        pScoreHat = self.ps.mean()
        return pScoreHat
        
    