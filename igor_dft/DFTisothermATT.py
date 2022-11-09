import math as mt
import numpy as np
import scipy as sp
from time import process_time
import matplotlib.pyplot as plt
from . import DFTpreprocessingATT as PP
from . import DFTcoreATT as CORE
from . import DFTinputsATT as IN
from . import DFTmathfunc as MATH

r = np.linspace(0+1E-3,PP.Hcc-1E-3,PP.NPmesh)
Vext = CORE.DFT_ExternalPotential(r)
rho0 = PP.rhob*np.ones(PP.NPmesh)


NPisotherm = 10
Press = np.linspace(0.001 , 0.060 , NPisotherm)

rhob   = np.zeros( NPisotherm )
mubres = np.zeros( NPisotherm )
gammaA = gammaD  = np.zeros( NPisotherm )
rhoA   = rhoD    = np.zeros((NPisotherm, PP.NPmesh))

print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXX ISOTHERM CONSTRUCTION XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX ASCENT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print(' ')
START_TIME = process_time()
for n in range(0,NPisotherm):
    print('=====================================================================================')
    print('==================================== POINT ',n,' =====================================')
    print('=====================================================================================')
    print(' ')
    [ rhob[n] , mubres[n] ] = PP.DFT_BulkEOS(PP.T,Press[n])
    if n == 0 :
        rho0 = rhob[0]*np.ones(PP.NPmesh)
    else :
        rho0 = rhoA[n-1,:]
    rhoA[n] = CORE.DFT_SolverPicard(r,rho0,PP.mubres[n],Vext,PP.alpha,PP.TOL,PP.ITmax)
    gammaA[n] = MATH.MATH_TrapezoidalRule(rhoA[n],r[0],r[-1])
    print(' ')
    print('=====================================================================================')
    print('=====================================================================================')
    print('=====================================================================================')
    print(' ')
END_TIME = process_time()
ASCENT_ELAPSED_TIME = END_TIME - START_TIME
print(' ')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('Elapsed time (h): ', ASCENT_ELAPSED_TIME/3600)
print(' ')

      
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXX ISOTHERM CONSTRUCTION XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX DESCENT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print(' ')
START_TIME = process_time()
for n in range(NPisotherm-1,-1,-1):
    print('=====================================================================================')
    print('==================================== POINT ',n,' =====================================')
    print('=====================================================================================')
    print(' ')
    if n == NPisotherm - 1 :
        rho0 = rhob[-1]*np.ones(PP.NPmesh)
    else :
        rho0 = rhoA[n+1]
    rhoA[n] = CORE.DFT_SolverPicard(r,rho0,PP.mubres[n],Vext,PP.alpha,PP.TOL,PP.ITmax)
    gammaA[n] = MATH.MATH_TrapezoidalRule(rhoA[n],r[0],r[-1])
    print(' ')
    print('=====================================================================================')
    print('=====================================================================================')
    print('=====================================================================================')
    print(' ')
END_TIME = process_time()
DESCENT_ELAPSED_TIME = END_TIME - START_TIME
print(' ')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
print('Elapsed time (h): ', DESCENT_ELAPSED_TIME/3600) 


