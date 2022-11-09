import math as mt
import numpy as np
import scipy.optimize as sp
from . import DFTinputsATT as IN

# CONSTANTE DE BOLTZMANN:
    
kB = IN.kB/IN.kB                 # (-)

# RELOADING INPUT VARIABLES:
    
P         = IN.P*(IN.sigmaff*1E-10)**3/(IN.epsilonff*IN.kB)
T         = IN.T/IN.epsilonff      
epsilonff = IN.epsilonff/IN.epsilonff
sigmaff   = IN.sigmaff/IN.sigmaff
epsilonsf = IN.epsilonsf/IN.epsilonff
sigmasf   = IN.sigmasf/IN.sigmaff
Hcc       = IN.Hcc/IN.sigmaff   
Delta     = IN.Delta/IN.sigmaff
rhos      = IN.rhos*IN.sigmaff**3

    
EXTERNAL_POTENTIAL_MODEL     = IN.EXTERNAL_POTENTIAL_MODEL
CHEMICAL_POTENTIAL_MODEL     = IN.CHEMICAL_POTENTIAL_MODEL
BULK_EOS                     = IN.BULK_EOS
FMT_FLAVOUR                  = IN.FMT_FLAVOUR
FLUID_ATTRACTIVE_INTERACTION = IN.FLUID_ATTRACTIVE_INTERACTION
ATTRACTIVE_POTENTIAL         = IN.ATTRACTIVE_POTENTIAL
    
alpha  = IN.alpha
TOL    = IN.TOL
ITmax  = IN.ITmax


# VARIABLES DEFINITIONS:

NPmesh = int((Hcc - 0)/IN.dz) + 1
EPS = 1E-2
zinf =   0 + EPS
zsup = Hcc - EPS

dhs = sigmaff
rmin = 2**(1/6)
rc   = 5.0

# DENSIDADE DO BULK:
   
def DFT_BulkEOS(T,P):
    kBT = kB*T
    if BULK_EOS == 'IG':
        rho = P/kBT
        mures = 0
    elif BULK_EOS == 'Carnahan-Starling':
        def EOSfobj(rho):
            eta = rho*(mt.pi*dhs**3)/6
            pHS = kBT*rho*(1 + eta + eta**2 - eta**3)/(1 - eta)**3
            RES = P - pHS
            return np.abs(RES)
        rho = sp.fsolve(EOSfobj, P/kBT)
        eta = mt.pi*rho*dhs**3/6
        mures = kBT*((8*eta - 9*eta**2 + 3*eta**3)/(1-eta)**3)
    elif BULK_EOS == 'BMCSL':
        def EOSfobj(rho):
            eta = rho*(mt.pi*dhs**3)/6
            pHS = kBT*rho*(1 + eta + eta**2 - eta**3)/(1 - eta)**3
            pWCA = 1/2*rho**2*( - 2*(1/2)*32/9*mt.pi*epsilonff*sigmaff**3 + 16/3*mt.pi*epsilonff*sigmaff**3*( (sigmaff/rc)**3 - 1/3*(sigmaff/rc)**9 ) )
            RES = P - (pHS + pWCA)
            return np.abs(RES)        
        rho = sp.fsolve(EOSfobj, P/kBT)
        eta = mt.pi*rho*dhs**3/6
        muHS = kBT*((8*eta - 9*eta**2 + 3*eta**3)/(1-eta)**3)
        mures = muHS  + rho*( - 2*(1/2)*32/9*mt.pi*epsilonff*sigmaff**3 + 16/3*mt.pi*epsilonff*sigmaff**3*( (sigmaff/rc)**3 - 1/3*(sigmaff/rc)**9 ) )
    return [rho , mures]


[rhob, mubres] = DFT_BulkEOS(T,P)   # (-)

###############################################################################