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

import sampleSimulation as ss
import numpy as np
import estimatorInterface as ei

# integration is numericially unstable, thus disable the warning!
import warnings
warnings.filterwarnings('ignore')

# load the data
Y,D,X,Z = ss.underiv_data(5000)

# Set the x for the MTE
x0 = np.mean(X, 0)
x0.shape = (1,2)

testEstObj = ei.treatmentEffectEst()

'''
**************************************************
Heck it: Linear Model, Normal Distribution
**************************************************
'''
# calculate the MTE
linear = 1
normal = 1
LLRMTE = 0
polyOrder = 4

testEstObj.assignData(Y, D, X, Z)
testEstObj.assignParam(linear, normal, LLRMTE, polyOrder)

testEstObj.estimate()


ateHat = testEstObj.ate(x0)
attHat = testEstObj.att(x0)
atutHat = testEstObj.atut(x0)

print "ATT , HUV algo 0.258; My algo ", "{0:.3f}".format(attHat)
print "ATUT, HUV algo 0.185; My algo ", "{0:.3f}".format(atutHat)
print "ATE , HUV algo 0.223; My algo ", "{0:.3f}".format(ateHat)


'''
**************************************************
Heckit on steroid: Linear Model, unknown Distribution
**************************************************
'''
# calculate the MTE
normal = 0
testEstObj.assignParam(linear, normal, LLRMTE, polyOrder)

testEstObj.estimate()

# Now set up the weight


ateHat = testEstObj.ate(x0)
attHat = testEstObj.att(x0)
atutHat = testEstObj.atut(x0)

print "ATT , HUV algo 0.279; My algo ", "{0:.3f}".format(attHat)
print "ATUT, HUV algo 0.128; My algo ", "{0:.3f}".format(atutHat)
print "ATE , HUV algo 0.202; My algo ", "{0:.3f}".format(ateHat)


'''
**************************************************
Nonlinear outcome equation, polynomial MTE estimator
**************************************************
'''
# calculate the MTE
linear = 0
testEstObj.assignParam(linear, normal, LLRMTE, polyOrder)

testEstObj.estimate()

ateHat = testEstObj.ate(x0)
attHat = testEstObj.att(x0)
atutHat = testEstObj.atut(x0)

print "ATT , HUV algo 0.261; My algo ", "{0:.3f}".format(attHat)
print "ATUT, HUV algo 0.158; My algo ", "{0:.3f}".format(atutHat)
print "ATE , HUV algo 0.208; My algo ", "{0:.3f}".format(ateHat)


'''
**************************************************
Nonlinear outcome equation, kernel regression MTE estimator
**************************************************
'''
# calculate the MTE
LLRMTE = 1
testEstObj.assignParam(linear, normal, LLRMTE, polyOrder)

testEstObj.estimate()

ateHat = testEstObj.ate(x0)
attHat = testEstObj.att(x0)
atutHat = testEstObj.atut(x0)

print "ATT , HUV algo 0.261; My algo ", "{0:.3f}".format(attHat)
print "ATUT, HUV algo 0.158; My algo ", "{0:.3f}".format(atutHat)
print "ATE , HUV algo 0.209; My algo ", "{0:.3f}".format(ateHat)