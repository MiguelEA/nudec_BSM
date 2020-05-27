# Last modified 27/05/2020. Miguel Escudero Abenza, miguel.escudero@kcl.ac.uk
# Check ArXiv:1812.05605 for the relevant equations
# Temperature evolution including a particle that is in thermal equilibrium with the neutrino sector
# The command to run should look like
# python WIMP_nu.py 12.9 2 F test_nu
# where 12.9 corresponds to the particle mass in MeV, 2 are the NDOF, and F corresponds to Fermi Statistics.
# This will output the temperature evolution to Results/WIMPS/test_nu.dat and Results/WIMPS/test_nu_Neff.dat the Neff value
#####################################################################
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad, odeint

# Flag to tell how to run the program
if len(sys.argv) < 5:
    print("ERROR"+"\n"+"NEED MORE INFORMATION, you should run for example:"+"\n"+"  python WIMP_nu.py 12.9 2 F test_nu")
    sys.exit()
else:
    MDM, gDM, SPIN, nameoffile  = float(sys.argv[1]), float(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4])

# WIMP Thermodynamics
if SPIN == 'B' or SPIN == 'BOSE' or SPIN == 'Bose':
    def rho_DM(T,MDM):
        if T < MDM/30.0: return 0.0
        else:            return gDM/(2*np.pi**2)*T**4*quad(lambda E: E**2*(E**2-(MDM/T)**2)**0.5/(np.exp(E)-1.) ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
    def p_DM(T,MDM):
        if T < MDM/30.0:  return 0.0
        else:           return gDM/(6*np.pi**2)*T**4*quad(lambda E: (E**2-(MDM/T)**2)**1.5/(np.exp(E)-1.) ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
    def drho_DMdT(T,MDM):
        if T < MDM/30.0:  return 0.0
        else:          return gDM/(2*np.pi**2)*T**3*quad(lambda E: 0.25*E**3*(E**2-(MDM/T)**2)**0.5*np.sinh(E/2.0)**-2 ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
if SPIN == 'F' or SPIN == 'Fermi' or SPIN == 'Fermion':
    def rho_DM(T,MDM):
        if T < MDM/30.0:            return 0.0
        else:            return gDM/(2*np.pi**2)*T**4*quad(lambda E: E**2*(E**2-(MDM/T)**2)**0.5/(np.exp(E)+1.) ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
    def p_DM(T,MDM):
        if T < MDM/30.0:            return 0.0
        else:            return gDM/(6*np.pi**2)*T**4*quad(lambda E: (E**2-(MDM/T)**2)**1.5/(np.exp(E)+1.) ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
    def drho_DMdT(T,MDM):
        if T < MDM/30.0:    return 0.0
        else:            return gDM/(2*np.pi**2)*T**3*quad(lambda E: 0.25*E**3*(E**2-(MDM/T)**2)**0.5*np.cosh(E/2.0)**-2 ,MDM/T,100,epsabs=1e-12,epsrel = 1e-12)[0]

#####################################################################
####################  RELEVANT CODE STARTS HERE  ####################
#####################################################################


# QED O(e^3)
P_int      =interp1d(np.loadtxt("QED/QED_P_int.cvs")[:,0],np.loadtxt("QED/QED_P_int.cvs")[:,1]+np.loadtxt("QED/QED_P_int.cvs")[:,2],bounds_error=False,fill_value=0.0,kind='linear')
dP_intdT   =interp1d(np.loadtxt("QED/QED_dP_intdT.cvs")[:,0],np.loadtxt("QED/QED_dP_intdT.cvs")[:,1]+np.loadtxt("QED/QED_dP_intdT.cvs")[:,2],bounds_error=False,fill_value=0.0,kind='linear')
d2P_intdT2 =interp1d(np.loadtxt("QED/QED_d2P_intdT2.cvs")[:,0],np.loadtxt("QED/QED_d2P_intdT2.cvs")[:,1]+np.loadtxt("QED/QED_d2P_intdT2.cvs")[:,2],bounds_error=False,fill_value=0.0,kind='linear')

## Uncomment in order to remove the QED corrections
#P_int, dP_intdT, d2P_intdT2   = lambda x: 0, lambda x: 0, lambda x: 0

# All in MeV Units!
GF  = 1.1663787e-5*1e-6 #in MeV^{-2}
me  = 0.511
Mpl = 1.22091e19*1e3

# Conversion factor to convert MeV^-1 into seconds
FAC = 1./(6.58212e-22)

# Left and Right nu-e couplings as relevant for E < 10 MeV. From the EW review of the PDG
geL, geR, gmuL, gmuR = 0.727, 0.233, -0.273, 0.233

# Thermodynamics
def rho_nu(T):  return 2 * 7./8. * np.pi**2/30. * T**4
def rho_gam(T): return 2 * np.pi**2/30. * T**4
def rho_e(T):
    if T < me/30.0:        return 0.0
    else:        return 4./(2*np.pi**2)*T**4*quad(lambda E: E**2*(E**2-(me/T)**2)**0.5/(np.exp(E)+1.) ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
def p_e(T):
    if T < me/30.0:        return 0.0
    else:        return 4./(6*np.pi**2)*T**4*quad(lambda E: (E**2-(me/T)**2)**1.5/(np.exp(E)+1.) ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]

# Derivatives
def drho_nudT(T):     return 4*rho_nu(T)/T
def drho_gamdT(T):    return 4*rho_gam(T)/T
def drho_edT(T):
    if T < me/30.0:  return 0.0
    else:        return 4./(2*np.pi**2)*T**3*quad(lambda E: 0.25*E**3*(E**2-(me/T)**2)**0.5*np.cosh(E/2.0)**-2 ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]

# Rho tot
def rho_tot(T_gam,T_nue,T_numu,MDM):
    return rho_gam(T_gam) + rho_e(T_gam) + rho_DM(T_nue,MDM) + rho_nu(T_nue) + 2*rho_nu(T_numu) - P_int(T_gam) + T_gam*dP_intdT(T_gam)
# P tot
def p_tot(T_gam,T_nue,T_numu,MDM):
    return  1./3. * rho_gam(T_gam) + p_e(T_gam) + p_DM(T_nue,MDM) + 1./3. * rho_nu(T_nue) + 1./3. * 2*rho_nu(T_numu) + P_int(T_gam)
# Hubble
def Hubble(T_gam,T_nue,T_numu,MDM):
    return FAC * (rho_tot(T_gam,T_nue,T_numu,MDM)*8*np.pi/(3*Mpl**2))**0.5
# Neff Definition
def Neff_func(T_gam,T_nue,T_numu):
    return 8./7.*(11./4.)**(4./3.)* (rho_nu(T_nue) + 2*rho_nu(T_numu)) / rho_gam(T_gam)

# Suppression of the rates as a result of a non-negligible electron mass
f_nue_s  = interp1d(np.loadtxt("SM_Rates/nue_scatt.dat")[:,0],np.loadtxt("SM_Rates/nue_scatt.dat")[:,1],kind='linear')
f_numu_s = interp1d(np.loadtxt("SM_Rates/numu_scatt.dat")[:,0],np.loadtxt("SM_Rates/numu_scatt.dat")[:,1],kind='linear')
f_nue_a  = interp1d(np.loadtxt("SM_Rates/nue_ann.dat")[:,0],np.loadtxt("SM_Rates/nue_ann.dat")[:,1],kind='linear')
f_numu_a = interp1d(np.loadtxt("SM_Rates/numu_ann.dat")[:,0],np.loadtxt("SM_Rates/numu_ann.dat")[:,1],kind='linear')

##Uncomment in order to remove the effect of m_e in the rates.
#f_nue_s, f_numu_s, f_nue_a, f_numu_a   = lambda T : 1, lambda T : 1, lambda T : 1, lambda T : 1

def Ffunc_nue_e(T1,T2):
    return 32* 0.884 *(T1**9 - T2**9) * f_nue_a(T1)  + 56 * 0.829 * f_nue_s(T1)  *T1**4*T2**4*(T1-T2)
def Ffunc_numu_e(T1,T2):
    return 32* 0.884 *(T1**9 - T2**9) * f_numu_a(T1) + 56 * 0.829 * f_numu_s(T1) *T1**4*T2**4*(T1-T2)
def Ffunc(T1,T2):
    return 32* 0.884 *(T1**9 - T2**9) + 56* 0.829 *T1**4*T2**4*(T1-T2)

# Energy Transfer Rates
def DeltaRho_nue(T_gam,T_nue,T_numu):
    return FAC * GF**2/np.pi**5 * ( 4* (geL**2 + geR**2) * Ffunc_nue_e(T_gam,T_nue)  + 2*Ffunc(T_numu,T_nue) )
def DeltaRho_numu(T_gam,T_nue,T_numu):
    return FAC * GF**2/np.pi**5 * ( 4* (gmuL**2 + gmuR**2) * Ffunc_numu_e(T_gam,T_numu) -   Ffunc(T_numu,T_nue) )

# Temperature Evolution Equations
def dTnu_dt(T_gam,T_nue,T_numu,MDM):
    return -(Hubble(T_gam,T_nue,T_nue,MDM) * (3 * 4 * rho_nu(T_nue) + 3*(rho_DM(T_nue,MDM)+p_DM(T_nue,MDM))) - (2*DeltaRho_numu(T_gam,T_nue,T_numu) + DeltaRho_nue(T_gam,T_nue,T_numu)))/ (3*drho_nudT(T_nue) + drho_DMdT(T_nue,MDM) )
def dTgam_dt(T_gam,T_nue,T_numu,MDM):
    return -(Hubble(T_gam,T_nue,T_nue,MDM)*( 4*rho_gam(T_gam) + 3*(rho_e(T_gam)+p_e(T_gam)) + 3 * T_gam * dP_intdT(T_gam)) + DeltaRho_nue(T_gam,T_nue,T_numu) + 2*DeltaRho_numu(T_gam,T_nue,T_numu) )/( drho_gamdT(T_gam) + drho_edT(T_gam)  + T_gam * d2P_intdT2(T_gam) )

def dT_totdt(vec,t,MDM):
    T_gam, T_nu = vec
    return [dTgam_dt(T_gam,T_nu,T_nu,MDM),dTnu_dt(T_gam,T_nu,T_nu,MDM)]


#Start the integration at a common temperature of 20 MeV, which corresponds to t ~ 2*10^{-3} s
T0 = 10.0
t0 = 1./(2*Hubble(T0,T0,T0,MDM))

# Finish the calculation at t = 5*10^4 seconds, which will correspond to T ~ 5*10^{-3} MeV
t_max = 5e4

# Calculate
tvec = np.logspace(np.log10(t0),np.log10(t_max),num=300)
sol = odeint(dT_totdt, [T0,T0], tvec, args=(MDM,), rtol = 1e-8, atol= 1e-8)

# Display Results
print("Neff      = ", round(Neff_func(sol[-1,0],sol[-1,1],sol[-1,1]),5))
print("Tgam/Tnu  = ", round(sol[-1,0]/sol[-1,1],5))

# Output Results
with open("Results/WIMPS/"+nameoffile+"_Neff.dat", 'w') as f:
    f.write(str(round(Neff_func(sol[-1,0],sol[-1,1],sol[-1,1]),5)))

# Calculate some Thermodynamic quantities
rho_vec, p_vec = np.zeros(len(tvec)),np.zeros(len(tvec))
for i in range(len(tvec)):
    rho_vec[i], p_vec[i]  = rho_tot(sol[i,0], sol[i,1], sol[i,1],MDM), p_tot(sol[i,0], sol[i,1], sol[i,1],MDM)

# Store results
print("Thermodynamics output to : Results/WIMPS/"+nameoffile+".dat")
with open("Results/WIMPS/"+nameoffile+".dat", 'w') as f:
    f.write("# Particle in thermal equilibrium with the neutrino sector of the plasma"+"\n")
    f.write("# Neff        = "),f.write("{:<12}".format(round(Neff_func(sol[-1,0],sol[-1,1],sol[-1,1]),5))),f.write("\n")
    f.write("# Tgam/Tnu    = "),f.write("{:<12}".format(round(sol[-1,0]/sol[-1,1],5))),f.write("\n")
    f.write("{:<13}".format("#T_gamma (MeV)")), f.write("{:<13}".format("T_gam/T_nu")),f.write("{:<14}".format("R_tot/T_gam^4 ")),f.write("{:<14}".format("P_tot/T_gam^4"+"\n"))
    for i in range(len(tvec)):
        f.write("{:.7E}".format(1e3*sol[i,0])),f.write("  "), f.write("{:<12}".format(round(sol[i,0]/sol[i,1],6))), f.write("{:.7E}".format(rho_vec[i]/sol[i,0]**4)),f.write("  "),f.write("{:.7E}".format(p_vec[i]/sol[i,0]**4)),f.write("\n")




