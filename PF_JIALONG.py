import numpy as np
import subprocess as sp
import sys 

def get_globals():
    global NF,NM,wc,A0,Epol
    NF = 5
    NM = 50
    wc = float (sys.argv[2]) / 27.2114 #eV
    A0 = float (sys.argv[1]) #0.1 #a.u.
    Epol = np.array([   1,0,0   ])
    #Epol = Epol / np.linalg.norm(Epol)

def get_data():
    Mu_data = np.loadtxt("../dipole_matrix_E.dat") 
    E_ad = np.loadtxt("../adiabtic_energies.dat") [:NM,1] / 27.2114 #eV to a.u.
    
    Mu_ad = np.zeros((NM, NM))
    for J,K,Dx,Dy,Dz in Mu_data:
        J = int (J)
        K = int (K)
        if (J <= NM-1):
            if (K <= NM-1):
                Mu_ad[J,K] = np.dot(Epol,[Dx,Dy,Dz])
                Mu_ad[K,J] = Mu_ad[J,K]


    return E_ad, Mu_ad

def build_H_PF(E_ad, Mu_ad):
    a_op = np.zeros((NF,NF))
    for n in range (1,NF):
        a_op[n,n-1]=np.sqrt(n)
    a_op = a_op.T

    I_m = np.identity(NM)
    I_ph = np.identity(NF)
    H_PF = np.kron(np.diag(E_ad), I_ph) #Hel
    H_PF += wc * np.kron(I_m, np.diag(np.arange(NF))) #+Hph
    H_PF += wc * A0 * np.kron(Mu_ad, a_op.T + a_op) #+Hel-ph
    H_PF += wc * A0**2 * np.kron(Mu_ad @ Mu_ad, I_ph) #+HDSE

    return H_PF

def save_results(E_PF, U_PF):
    sp.call("mkdir -p data_PF", shell=True)
    np.savetxt(f"data_PF/Epol_A0{round(A0,4)}_wc{round(wc*27.2114,4)}_NM{NM}_NF{NF}.dat",E_PF*27.2114)
    np.savetxt(f"data_PF/Upol_A0{round(A0,4)}_wc{round(wc*27.2114,4)}_NM{NM}_NF{NF}.dat",U_PF)
    #np.save(f"data_PF/Upol_A0{round(A0,4)}_wc{round(wc*27.2114,4)}_NM{NM}_NF{NF}.dat.npy",U_PF)

def main():
    get_globals()
    E_ad, Mu_ad = get_data() #Step1
    H_PF = build_H_PF(E_ad, Mu_ad) #Step2
    E_PF, U_PF = np.linalg.eigh(H_PF) #Step3
    #print (E_PF.shape, U_PF.shape)
    save_results(E_PF, U_PF) #Step4

if (__name__ ==  "__main__"):
    main()