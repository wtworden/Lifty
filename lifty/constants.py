
from lifty.sage_types import *

z = polygen(QQ)

roots = complex_roots((z**2+z)**2+z, retval='algebraic')

RABBIT = roots[3][0]

CORABBIT = roots[2][0]

AIRPLANE = roots[0][0]