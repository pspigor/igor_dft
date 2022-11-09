# CONSTANTE DE BOLTZMANN:

kB = 1.380649E-23

# CONDIÇÕES NO BULK DO SISTEMA:

P = 0.7*1.01325E+5         # (Pa)
T = 77.4                   # (K)


# PARÂMETROS DO POTENCIAL INTERMOLECULAR:

sigmaff   = 3.575          # (Å)
epsilonff = 94.45          # (K)
sigmasf   = 3.494          # (Å)
epsilonsf = 53.22          # (K)
rhos      = 0.114          # (1/Å³)

# DISTÂNCIA ENTRE AS PAREDES:

Hcc   = 36.00              # (Å) 
Hcc   = 360.00
Delta =  3.35              # (Å)

# SELEÇÃO DO MODELO:
    
EXTERNAL_POTENTIAL_MODEL     = 'Steele'
CHEMICAL_POTENTIAL_MODEL     = 'BMCSL'
BULK_EOS                     = 'BMCSL'
FMT_FLAVOUR                  = 'WB'
FLUID_ATTRACTIVE_INTERACTION = 'MeanField'
ATTRACTIVE_POTENTIAL         = 'WCA'

# PARÂMETROS DE CONTROLE COMPUTACIONAL:
    
alpha = 0.0005
TOL   = 1E-6
ITmax = 300
#ITmax = 10000
dz    = 1E-2

###############################################################################
