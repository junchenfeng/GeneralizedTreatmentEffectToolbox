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
This script estimates the treatment effect weights, according to the formula listed on table 1B page 396
available weights include
ate = 1
tt = \int_u^1 f(p|x) dp * (1/E(P|x))
tut = \int_0^u f(p|x) dp * (1/E(1-P|x))

prte: to come
IV: to come
'''


#from scipy import integrate

class weightGenerator():

    def __init__(self):
        self.pScoreObj =[]
        
    def _assignPScoreFnc(self, pScoreObj):
        self.pScoreObj = pScoreObj
        
    def wATE(self, x, u):
        weight = 1
        return weight
    
    def wTT(self, x, u):
        # anti-bug
        assert(isinstance(u, float))
        assert(u>=0 and u<=1)
        
        # set up the function 
        #psIntFnc = lambda U: self.pScoreObj.calcfPs(x, U)
        # integrate 
        #weight = integrate.quad(psIntFnc, u ,1)[0] / self.pScoreObj.calcPs(x)

        '''
        The current method does not actually conditions on X! 
        Especially, it does not handle the fact that X could overlap with Z
        '''
        weight =  self.pScoreObj.calcfPs(x, u) / self.pScoreObj.calcPs(x)
        
        return weight
      
    def wTUT(self,x, u):
        # anti-bug
        assert(isinstance(u, float))
        assert(u>=0 and u<=1)
        
        # set up the function 
        #psIntFnc = lambda U: self.pScoreObj.calcfPs(x,U)   
        
        # integrate 
        #weight = integrate.quad(psIntFnc, 0 ,u)[0]  / (1 - self.pScoreObj.calcPs(x))
        weight =  (1-self.pScoreObj.calcfPs(x, u)) / (1-self.pScoreObj.calcPs(x))
        return weight

''' 
experiment with integration
'''
