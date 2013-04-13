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
import kernel as kw
'''
check kernelInt 
'''

assert kw.kernelInt(0) == 1
assert kw.kernelInt(-0.9) == 1
assert kw.kernelInt(0.9) == 1
assert kw.kernelInt(1.1) == 0
assert kw.kernelInt(-1.1) == 0
# test for longer arrays
x = np.array([0,-1.1,0.8])
trueVal = np.array([1,0,1], dtype = int)
z = kw.kernelInt(x)
assert (z == trueVal).any

'''
check tricube
'''
x = np.array([-1.1, -1., -0.5, 0, 0.5, 1, 1.1] )
z = kw.tricube(x, 0, 1)
trueVal = np.array([0, 0, (7./8.)**3*(70./81.), (70./81.),(7./8.)**3*(70./81.),0,0])
assert (z == trueVal).any

print 'Kernel passes unit test.'
