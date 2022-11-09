import math as mt
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from . import DFTpreprocessingATT as PP
from . import DFTcoreATT as CORE
from . import DFTinputsATT as IN


r = np.linspace(0+1E-3,PP.Hcc-1E-3,PP.NPmesh)
Vext = CORE.DFT_ExternalPotential(r)
rho0 = PP.rhob*np.ones(PP.NPmesh)
rho = CORE.DFT_SolverPicard(r,rho0,PP.mubres,Vext,PP.alpha,PP.TOL,PP.ITmax)
dAres_drho = CORE.DFT_ResidualHelmohotzFunctionalDerivative(r,rho)

# plt.plot(r,Vext,color = 'orange', linewidth = 0.9)
# plt.plot(r,dAres_drho,color = 'red', linewidth = 0.9)
# plt.plot(r,PP.mubres*np.ones(PP.NPmesh),color = 'blue', linewidth = 0.9)
# plt.legend(['Potencial externo','Funcional (NLDFT-FMT)','Potencial químico'], fontsize = 14)   
# plt.title("Picard terms comparison", fontsize = 15)
# plt.xlabel('$ z \; (Å)$', fontsize = 14)
# plt.ylabel('$V^{\mathrm{ext}}$', fontsize = 14)
# plt.tight_layout()
# plt.show()

MC_DATA = np.loadtxt('MC_N2.dat')

#plt.plot(MC_DATA[:,0]/IN.sigmaff, MC_DATA[:,1]*IN.sigmaff**3, 'ro', markersize=7, markeredgewidth=1.3, mfc = 'none')
plt.plot(r, rho*(PP.sigmaff)**3,color = 'blue', linewidth = 1.0)
#plt.legend(['GCMC','NLDFT-FMT'], fontsize = 14)   
#plt.title("Density distribution", fontsize = 15)
plt.xlabel('$ z/ \sigma$', fontsize = 14)
plt.ylabel('$ρ \, \sigma^3 $', fontsize = 14)
plt.xlim((0.0,PP.Hcc/2))
plt.ylim((0.0,9.0))
plt.tight_layout()
plt.show()  