import numpy as np
from glob import glob

dA = 0.2
THETA_LIST = np.arange( 0, np.pi+dA, dA )
PHI_LIST   = np.arange( 0, np.pi+dA, dA )

A0 = 0.2
WC = 1.5

##### optimal theta/phi #####
E = np.zeros( (2, len(THETA_LIST), len(PHI_LIST)) ) 
for THETAi,THETA in enumerate( THETA_LIST ):
    for PHIi,PHI in enumerate( PHI_LIST ):
        THETA = round(THETA,2)
        PHI   = round(PHI,2)
        E[0,THETAi, PHIi] = np.loadtxt( f"TD_R/PF_NF10_NM50_THETA_PHI/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_10_NM_50.dat" )[0] / 27.2114 * 630
        E[1,THETAi, PHIi] = np.loadtxt( f"TD_TS/PF_NF10_NM50_THETA_PHI/data_PF/E_THETA_{THETA}_PHI_{PHI}_A0_{A0}_WC_{WC}_NF_10_NM_50.dat" )[0] / 27.2114 * 630

dE = E[1,:,:] - E[0,:,:]
CRITICAL_ANGLE_MIN = np.unravel_index( np.argmin(dE), dE.shape )
CRITICAL_ANGLE_MAX = np.unravel_index( np.argmax(dE), dE.shape )

print("Critical Point MIN:")
print( "\tTheta = %1.1f (deg.)" % (THETA_LIST[CRITICAL_ANGLE_MIN[0]] * 180/np.pi) )
print( "\tPhi   = %1.1f (deg.)" % (PHI_LIST[CRITICAL_ANGLE_MIN[1]] * 180/np.pi) )

print("Critical Point MAX:")
print( "\tTheta = %1.1f (deg.)" % (THETA_LIST[CRITICAL_ANGLE_MAX[0]] * 180/np.pi) )
print( "\tPhi   = %1.1f (deg.)" % (PHI_LIST[CRITICAL_ANGLE_MAX[1]] * 180/np.pi) )
