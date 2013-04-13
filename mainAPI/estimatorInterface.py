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

from scipy import integrate
import mteObj as mo
import weightObj as wo

class treatmentEffectEst():
    
    def __init__(self):
        self.mteObj = mo.hackit()
        self.weightObj = wo.weightGenerator()
    
    
    def assignData(self,Y, D, X, Z):
        
        self.mteObj.assignData(Y, D, X, Z)
        
        
    def assignParam(self,linear, normal, LLRMTE, polyOrder):
        self.mteObj.assignParam(linear, normal, LLRMTE, polyOrder)
        
        
    def estimate(self):
        self.mteObj.estimateParam()        
        self.weightObj._assignPScoreFnc(self.mteObj)
  

    def ate(self, x):
        deltaHat = integrate.quad(lambda U: self.weightObj.wATE(x, U) * self.mteObj.mtePredict(x, U), 0, 1)
        return deltaHat[0]
    
    def att(self, x):
        deltaHat = integrate.quad(lambda U: self.weightObj.wTT(x, U) * self.mteObj.mtePredict(x, U), 0, 1)
        return deltaHat[0]  
    
    def atut(self, x):
        deltaHat = integrate.quad(lambda U: self.weightObj.wTUT(x, U) * self.mteObj.mtePredict(x, U), 0, 1)
        return deltaHat[0]   