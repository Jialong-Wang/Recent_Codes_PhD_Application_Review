import numpy as np
from matplotlib import pyplot as plt
 
NF = 5
NM = 50
Pol = "y"

A0_list = np.arange(0.0,0.1+0.01,0.01)
WC_list = np.arange(0.0,10+1,1)

Na0 = len (A0_list)
NWC = len (WC_list)
N_pol = NF * NM
E_pol = np.zeros ((3,N_pol,Na0,NWC))

for A0_ind, A0 in enumerate(A0_list):
    for WC_ind, WC in enumerate(WC_list):
        A0 = round(A0, 4)
        WC = round(WC, 4)
        E_pol [0,:,A0_ind,WC_ind] = 1/27.2114*630*np.loadtxt(f"../3_TD_R/PF_NM{NM}_NF{NF}_E{Pol}/data_PF/Epol_A0{round(A0,4)}_wc{round(WC,4)}_NM{NM}_NF{NF}.dat")
        E_pol [1,:,A0_ind,WC_ind] = 1/27.2114*630*np.loadtxt(f"../3_TD_TS/PF_NM{NM}_NF{NF}_E{Pol}/data_PF/Epol_A0{round(A0,4)}_wc{round(WC,4)}_NM{NM}_NF{NF}.dat")
        E_pol [2,:,A0_ind,WC_ind] = 1/27.2114*630*np.loadtxt(f"../3_TD_P/PF_NM{NM}_NF{NF}_E{Pol}/data_PF/Epol_A0{round(A0,4)}_wc{round(WC,4)}_NM{NM}_NF{NF}.dat")

        #print (A0,WC, E_pol[1,0,A0_ind,WC_ind] - E_pol[0,0,A0_ind,0], E_pol [2,0,A0_ind,WC_ind] - E_pol [0,0,A0_ind,0] )
        print ("A0,WC:",A0,WC)

print (len(A0_list),len(WC_list))
print (Na0,NWC)

#plot polariton energy as a function of cavity coupling
E_0 =  E_pol [0,0,0,0]

plt.plot(A0_list,E_pol[0,1,:,49*NWC//100] - E_0)
plt.plot(A0_list,E_pol[0,2,:,49*NWC//100] - E_0)
plt.plot(A0_list,E_pol[0,3,:,49*NWC//100] - E_0)
plt.plot(A0_list,E_pol[0,4,:,49*NWC//100] - E_0)
plt.plot(A0_list,E_pol[0,5,:,49*NWC//100] - E_0)
plt.savefig(f"E_A0_E{Pol}.jpg")
plt.clf()

#plot polariton energy as a function of cavity coupling
E_0 =  E_pol [0,0,0,0]

print (E_pol.shape)
plt.plot(WC_list,E_pol[0,1,Na0-1,:] - E_0)
plt.plot(WC_list,E_pol[0,2,Na0-1,:] - E_0)
plt.plot(WC_list,E_pol[0,3,Na0-1,:] - E_0)
plt.plot(WC_list,E_pol[0,4,Na0-1,:] - E_0)
plt.plot(WC_list,E_pol[0,5,Na0-1,:] - E_0)
plt.savefig(f"E_WC_E{Pol}.jpg")
plt.clf()

#plot ground state polariton energy for TS as the function of coupling strength for different wc

for WC_ind, WC in enumerate(WC_list):
    WC = round(WC, 4)
    E_diff = E_pol [1,0,:,WC_ind] - E_pol [0,0,:,WC_ind]
    plt.plot(A0_list,E_diff, label= f"$\omega_c$ = {WC} eV")

plt.legend()
plt.savefig(f"E_TS-R_A0_E{Pol}.jpg")
plt.clf()

#plot ground state polariton energy for product as the function of coupling strength for different wc

for WC_ind, WC in enumerate(WC_list):
    WC = round(WC, 4)
    E_diff = E_pol [2,0,:,WC_ind] - E_pol [0,0,:,WC_ind]
    plt.plot(A0_list,E_diff, label= f"$\omega_c$ = {WC} eV")

plt.legend()
plt.savefig(f"E_P-R_A0_E{Pol}.jpg")

#Maximum TS shift
print ("Maximum TS shift")
print ("A0",A0_list[0],"WC",WC_list[-1],"0",E_pol [1,0,0,-1] - E_pol [0,0,0,-1])
print ("A0",A0_list[-1],"WC",WC_list[-1],"0",E_pol [1,0,-1,-1] - E_pol [0,0,-1,-1])

#Maximum P shift
print ("Maximum P shift")
print ("A0",A0_list[0],"WC",WC_list[-1],"0",E_pol [2,0,0,-1] - E_pol [0,0,0,-1])
print ("A0",A0_list[-1],"WC",WC_list[-1],"0",E_pol [2,0,-1,-1] - E_pol [0,0,-1,-1])