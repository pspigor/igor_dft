import math as mt
import numpy as np
import scipy as sp
from time import process_time
from . import DFTpreprocessingATT as PP
from . import DFTmathfunc as MT
import matplotlib.pyplot as plt


def DFT_LocalDensity(r,rho,mubres,Vext):
    kBT = PP.kB*PP.T
    #mub_res    = DFT_ResidualChemicalPotential(PP.T,PP.P,rho)
    dAres_drho = DFT_ResidualHelmohotzFunctionalDerivative(r,rho)
    rho = PP.rhob*np.exp(1/kBT*( mubres - dAres_drho - Vext )) 
    return rho 

def DFT_SolverPicard(r,rho0,mubres,Vext,alpha,TOL,ITmax):
    COND = True
    ERRO = np.zeros(PP.NPmesh)
    IT = 1
    START_TIME = process_time()
    while COND == True:
        if IT in range(1,50):
            Alpha = alpha
        else :
            Alpha = 50*alpha
        rho_calc = DFT_LocalDensity(r,rho0,mubres,Vext)
        rho = (1 - Alpha)*rho0 + Alpha*rho_calc
        ERRO = np.abs(rho - rho0)
        maxERRO = np.max(ERRO)
        sumERRO = np.sum(ERRO)
        COND = sumERRO > TOL and IT < ITmax
        rho0 = np.copy(rho)
        print('###################################### PICARD #########################################')
        print('IT: ',IT)
        print('α: ', Alpha)
        print('RES: ', sumERRO)
        if IT >= ITmax :
            print('Error: Maximum Number of Iterations Exceeded') 
        IT = IT + 1
        END_TIME = process_time()
    print('#######################################################################################')
    print('Elapsed time (min): ', (END_TIME-START_TIME)/60)    
    return rho

def DFT_ExternalPotential(r):
    if   PP.EXTERNAL_POTENTIAL_MODEL == 'NONE':
        V1 = V2 = np.zeros(PP.NPmesh)
        Vext = V1 + V2
        Vext[:50] = Vext[PP.NPmesh-50:] = 100
    elif PP.EXTERNAL_POTENTIAL_MODEL == 'Steele':
        z = r
        V1 = 2*mt.pi*PP.rhos*PP.epsilonsf*PP.sigmasf**2*PP.Delta*(2/5*(PP.sigmasf/z)**10 - (PP.sigmasf/z)**4 - 1/3*PP.sigmasf**4/(PP.Delta*(0.61*PP.Delta+z)**3) ) 
        z = PP.Hcc - r
        V2 = 2*mt.pi*PP.rhos*PP.epsilonsf*PP.sigmasf**2*PP.Delta*(2/5*(PP.sigmasf/z)**10 - (PP.sigmasf/z)**4 - 1/3*PP.sigmasf**4/(PP.Delta*(0.61*PP.Delta+z)**3) ) 
        Vext = V1 + V2
        Vext[Vext > 1E+2] = 1E+2
    return Vext

def DFT_FMTweight(z):
    R = 0.5*PP.dhs
    omega2 = 2*mt.pi*R*np.heaviside(R - np.abs(z),1)
    omega0 = omega2/(4*mt.pi*R**2)
    omega1 = omega2/(4*mt.pi*R)
    omega3 = mt.pi*(R**2 - z**2)*np.heaviside(R - np.abs(z),1)
    OMEGA2 = 2*mt.pi*z*np.heaviside(R - np.abs(z),1)
    OMEGA1 = OMEGA2/(4*mt.pi*R)
    return [ omega0, omega1, omega2 , omega3, OMEGA1 , OMEGA2 ]
    

def DFT_FMTpotential(r,rho):
    z = np.copy(r)
    NP = np.size(rho)
    omega0 = np.zeros(NP)
    omega1 = np.zeros(NP)
    omega2 = np.zeros(NP)
    omega3 = np.zeros(NP)
    OMEGA1 = np.zeros(NP)
    OMEGA2 = np.zeros(NP)
    n0     = np.zeros(NP)
    n1     = np.zeros(NP)
    n2     = np.zeros(NP)
    n3     = np.zeros(NP)
    N1     = np.zeros(NP)
    N2     = np.zeros(NP)
    for n in range(0,NP):
        [ omega0, omega1, omega2 , omega3 , OMEGA1 , OMEGA2 ] = DFT_FMTweight(z[n]*np.ones(NP) - z)
        n0[n] = MT.MATH_TrapezoidalRule(rho*omega0, PP.zinf, PP.zsup)
        n1[n] = MT.MATH_TrapezoidalRule(rho*omega1, PP.zinf, PP.zsup)
        n2[n] = MT.MATH_TrapezoidalRule(rho*omega2, PP.zinf, PP.zsup)
        n3[n] = MT.MATH_TrapezoidalRule(rho*omega3, PP.zinf, PP.zsup)
        N1[n] = MT.MATH_TrapezoidalRule(rho*OMEGA1, PP.zinf, PP.zsup)
        N2[n] = MT.MATH_TrapezoidalRule(rho*OMEGA2, PP.zinf, PP.zsup)
    if PP.FMT_FLAVOUR == 'WB' :
        PHI = - n0*np.log(1 - n3) + (n1*n2 - N1*N2)/(1 - n3) + (n2**3 - 3*n2*N2*N2)*(n3 + (1- n3)**2*np.log(1 - n3))/(36*np.pi*n3**2*(1 - n3)**2)
        dPHI_dn0 = - np.log(1 - n3)
        dPHI_dn1 = n2/(1 - n3)
        dPHI_dn2 = n1/(1 - n3) + 1/36*((3*n2**2 - 3*N2*N2)*(n3 + (1 - n3)**2*np.log(1 - n3)))/(np.pi*n3**2*(1 - n3)**2)
        dPHI_dn3 = n0/(1 - n3) + (n1*n2 - N1*N2)/(1 - n3)**2 - (n2**3 - 3*n2*N2*N2)/(36*np.pi*n3**3*(1 - n3)**3)*( n3*(n3**2 - 5*n3 + 2) + 2*(1 - n3)**3*np.log(1 - n3) )
        dPHI_dN1 = - N2/(1 - n3)
        dPHI_dN2 = - N1/(1 - n3) - 1/6*(n2*N2*(n3 + (1 - n3)**2*np.log(1 - n3)))/(np.pi*n3**2*(1 - n3)**2)
    elif PP.FMT_FLAVOUR == 'RF' :
        PHI = - n0*np.log(1 - n3) + (n1*n2  - N1*N2)/(1 - n3) + (n2**3  - 3*n2*N2*N2)/(24*np.pi*(1 - n3)**2)
        dPHI_dn0 = - np.log(1 - n3)
        dPHI_dn1 = n2/(1 - n3)
        dPHI_dn2 = n1/(1 - n3) + 1/24*(3*n2**2 - 3*N2*N2)/(np.pi*(1 - n3)**2)
        dPHI_dn3 = N1/(1 - n3)+ (n1*n2 - N1*N2)/(1 - n3)**2 + 1/12*(n2**3 - 3*n2*N2*N2)/(np.pi*(1 - n3)**3)
        dPHI_dN1 = - N2/(1 - n3)
        dPHI_dN2 = - N1/(1 - n3) - 1/4*(n2*N2)/(np.pi*(1 - n3)**2)
    return [ PHI , dPHI_dn0 , dPHI_dn1 , dPHI_dn2 , dPHI_dn3 , dPHI_dN1 , dPHI_dN2 ]

def DFT_HardSpheresHelmohotzFunctional(r,rho):
    kBT = PP.kB*PP.T
    [ PHI , dPHI_dn0 , dPHI_dn1 , dPHI_dn2 , dPHI_dn3 , dPHI_dN1 , dPHI_dN2 ] = DFT_FMTpotential(r,rho)
    Ahs = kBT*MT.MATH_TrapezoidalRule(PHI, PP.zinf, PP.zsup)
    return Ahs

def DFT_HardSpheresHelmohotzFunctionalDerivative(r,rho):
    kBT = PP.kB*PP.T
    z = np.copy(r)
    NP = np.size(rho)
    INTdPHI_dn0omega0 = np.zeros(NP)
    INTdPHI_dn1omega1 = np.zeros(NP)
    INTdPHI_dn2omega2 = np.zeros(NP)
    INTdPHI_dn3omega3 = np.zeros(NP)
    INTdPHI_dN1OMEGA1 = np.zeros(NP)
    INTdPHI_dN2OMEGA2 = np.zeros(NP)
    [ PHI , dPHI_dn0 , dPHI_dn1 , dPHI_dn2 , dPHI_dn3 , dPHI_dN1 , dPHI_dN2 ] = DFT_FMTpotential(r,rho)
    for n in range(0,NP):
        [ omega0, omega1, omega2 , omega3 , OMEGA1 , OMEGA2 ] = DFT_FMTweight(z[n]*np.ones(NP) - z)
        INTdPHI_dn0omega0[n] = MT.MATH_TrapezoidalRule(dPHI_dn0*omega0, PP.zinf, PP.zsup)
        INTdPHI_dn1omega1[n] = MT.MATH_TrapezoidalRule(dPHI_dn1*omega1, PP.zinf, PP.zsup)
        INTdPHI_dn2omega2[n] = MT.MATH_TrapezoidalRule(dPHI_dn2*omega2, PP.zinf, PP.zsup)
        INTdPHI_dn3omega3[n] = MT.MATH_TrapezoidalRule(dPHI_dn3*omega3, PP.zinf, PP.zsup)
        INTdPHI_dN1OMEGA1[n] = MT.MATH_TrapezoidalRule(dPHI_dN1*OMEGA1, PP.zinf, PP.zsup)
        INTdPHI_dN2OMEGA2[n] = MT.MATH_TrapezoidalRule(dPHI_dN2*OMEGA2, PP.zinf, PP.zsup)
    dAhs_drho = kBT*(INTdPHI_dn0omega0 + INTdPHI_dn1omega1 + INTdPHI_dn2omega2 + INTdPHI_dn3omega3 - INTdPHI_dN1OMEGA1 - INTdPHI_dN2OMEGA2)
    return dAhs_drho

def DFT_AttractivePotential(z):
    if PP.ATTRACTIVE_POTENTIAL == 'WCA' :
        NP = np.size(z)
        XI = np.zeros(NP)
        XI[np.abs(z) >= PP.rc  ] = 0
        XI[np.abs(z) >= PP.rmin] = 2/5*PP.epsilonff/z[np.abs(z) >= PP.rmin]**10 - PP.epsilonff/z[np.abs(z) >= PP.rmin]**4 - 2/5*PP.epsilonff/PP.rc**10 + PP.epsilonff/PP.rc**4
        XI[np.abs(z) <  PP.rmin] = 1/2*PP.epsilonff*(np.abs(z[np.abs(z) <  PP.rmin])**2 - PP.rmin**2) + 2/5*PP.epsilonff/PP.rmin**10 - PP.epsilonff/PP.rmin**4 - 2/5*PP.epsilonff/PP.rc**10 + PP.epsilonff/PP.rc**4
    return XI

def DFT_AttractiveHelmohotzFunctionalDerivative(r,rho):
    z = np.copy(r)
    NP = np.size(rho)
    dAatt_drho = np.zeros(NP)
    if PP.FLUID_ATTRACTIVE_INTERACTION == 'MeanField':
        for n in range(0,NP):
            XI = DFT_AttractivePotential(z[n]*np.ones(NP) - z)
            dAatt_drho[n] = MT.MATH_TrapezoidalRule(2*mt.pi*rho*XI, PP.zinf, PP.zsup)
    return dAatt_drho
    

def DFT_ResidualHelmohotzFunctionalDerivative(r,rho):
    if PP.BULK_EOS == 'IG':
        dAres_drho = 0
    elif PP.BULK_EOS == 'Carnahan-Starling':
        dAhs_drho = DFT_HardSpheresHelmohotzFunctionalDerivative(r,rho)
        dAres_drho = dAhs_drho
    elif PP.BULK_EOS == 'BMCSL':
        dAatt_drho = DFT_AttractiveHelmohotzFunctionalDerivative(r,rho)
        dAhs_drho = DFT_HardSpheresHelmohotzFunctionalDerivative(r,rho)
        dAres_drho = dAatt_drho + dAhs_drho
    return dAres_drho

def TESTING_IdealGasDensityProfile(r):
    kBT = PP.kB*PP.T
    Vext = DFT_ExternalPotential(r)
    rho = PP.rhob*np.exp(-Vext/kBT)
    return rho
    

