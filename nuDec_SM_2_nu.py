# Last modified 13/12/18. Miguel Escudero Abenza, miguel.escudero@kcl.ac.uk 
# Check ArXiv:1812.05605 for the relevant equations
# Calculation of the temperature evolution in the SM and Neff in the MB approximation
# Evolving separately the temperatures of nu_e and nu_mu
# Outputs a file with the temperature evolution as a function of time, named "SM/SM_emu.dat"
#####################################################################
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad, odeint

# QED O(alpha)
P_int      =interp1d(np.loadtxt("QED/QED_P_int.cvs")[:,0],np.loadtxt("QED/QED_P_int.cvs")[:,1],bounds_error=False,fill_value=0.0,kind='linear')
dP_intdT   =interp1d(np.loadtxt("QED/QED_dP_intdT.cvs")[:,0],np.loadtxt("QED/QED_dP_intdT.cvs")[:,1],bounds_error=False,fill_value=0.0,kind='linear')
d2P_intdT2 =interp1d(np.loadtxt("QED/QED_d2P_intdT2.cvs")[:,0],np.loadtxt("QED/QED_d2P_intdT2.cvs")[:,1],bounds_error=False,fill_value=0.0,kind='linear')

## Uncomment in order to remove the QED corrections
#P_int, dP_intdT, d2P_intdT2   = lambda x: 0, lambda x: 0, lambda x: 0

# All in GeV Units!
GF  = 1.1663787e-5 #in GeV^{-2}
me  = 0.511e-3
Mpl = 1.22091e19

# Conversion factor to convert GeV^-1 to seconds
FAC = 1./(6.58e-25)

# SM, neutrino couplings sW2 =  1-Mw^2/Mz^2
sW2 = 0.223
CAe, CVe, CAmu, CVmu = 0.5, 0.5 + 2*sW2, -0.5, -0.5 + 2*sW2

# Thermodynamics
def rho_nu(T):  return 2 * 7./8. * np.pi**2/30. * T**4
def rho_gam(T):  return 2 * np.pi**2/30. * T**4
def rho_e(T):
    if T < me/30.0: return 0.0
    else:        return 4./(2*np.pi**2)*T**4*quad(lambda E: E**2*(E**2-(me/T)**2)**0.5/(np.exp(E)+1.) ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]
def p_e(T):
    if T < me/30.0:        return 0.0
    else:        return 4./(6*np.pi**2)*T**4*quad(lambda E: (E**2-(me/T)**2)**1.5/(np.exp(E)+1.) ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]

# Derivatives
def drho_nudT(T):    return 4*rho_nu(T)/T
def drho_gamdT(T):    return 4*rho_gam(T)/T
def drho_edT(T):
    if T < me/30.0:        return 0.0
    else:        return 4./(2*np.pi**2)*T**3*quad(lambda E: 0.25*E**3*(E**2-(me/T)**2)**0.5*np.cosh(E/2.0)**-2 ,me/T,100,epsabs=1e-12,epsrel = 1e-12)[0]


# Rho tot
def rho_tot(T_gam,T_nue,T_numu):
    return rho_gam(T_gam) + rho_e(T_gam) + rho_nu(T_nue) + 2*rho_nu(T_numu) - P_int(T_gam) + T_gam*dP_intdT(T_gam)
# P tot
def p_tot(T_gam,T_nue,T_numu):
    return  1./3. * rho_gam(T_gam) + p_e(T_gam) + 1./3. * rho_nu(T_nue) + 1./3. * 2*rho_nu(T_numu) + P_int(T_gam)
# Hubble
def Hubble(T_gam,T_nue,T_numu):
    return FAC * (rho_tot(T_gam,T_nue,T_numu)*8*np.pi/(3*Mpl**2))**0.5

# Neff Definition
def Neff_func(T_gam,T_nue,T_numu):
    return 8./7.*(11./4.)**(4./3.)* (rho_nu(T_nue) + 2*rho_nu(T_numu)) / rho_gam(T_gam)

# Energy Transfer Rates
def DeltaRho_nue(T_gam,T_nue,T_numu):
    return FAC * 2 * GF**2*(CAe**2 + CVe**2)/np.pi**5 * (16*(T_gam**9 - T_nue**9) + 7 * T_gam**4 * T_nue**4 * (T_gam-T_nue) ) + FAC  * GF**2/np.pi**5 * (16*(T_numu**9 - T_nue**9) + 7 * T_numu**4 * T_nue**4 * (T_numu-T_nue) )
def DeltaRho_numu(T_gam,T_nue,T_numu):
    return FAC * 2 * GF**2*(CAmu**2 + CVmu**2)/np.pi**5 * (16*(T_gam**9 - T_numu**9) + 7 * T_gam**4 * T_numu**4 * (T_gam-T_numu) ) - FAC * 0.5 * GF**2/np.pi**5 * (16*(T_numu**9 - T_nue**9) + 7 * T_numu**4 * T_nue**4 * (T_numu-T_nue) )

# Temperature Evolution Equations
def dTnue_dt(T_gam,T_nue,T_numu):
    return -4*Hubble(T_gam,T_nue,T_numu) * rho_nu(T_nue) / drho_nudT(T_nue) + DeltaRho_nue(T_gam,T_nue,T_numu)/ drho_nudT(T_nue)

def dTnumu_dt(T_gam,T_nue,T_numu):
    return -4*Hubble(T_gam,T_nue,T_numu) * rho_nu(T_numu) / drho_nudT(T_numu) + 2*DeltaRho_numu(T_gam,T_nue,T_numu)/(2* drho_nudT(T_numu))

def dTgam_dt(T_gam,T_nue,T_numu):
    return -(Hubble(T_gam,T_nue,T_numu)*( 4*rho_gam(T_gam) + 3*(rho_e(T_gam)+p_e(T_gam)) + 3 * T_gam * dP_intdT(T_gam)) + DeltaRho_nue(T_gam,T_nue,T_numu) + 2*DeltaRho_numu(T_gam,T_nue,T_numu) )/( drho_gamdT(T_gam) + drho_edT(T_gam) + T_gam * d2P_intdT2(T_gam) )

def dT_totdt(vec,t):
    T_gam, T_nue, T_numu = vec
    return [dTgam_dt(T_gam,T_nue,T_numu),dTnue_dt(T_gam,T_nue,T_numu),dTnumu_dt(T_gam,T_nue,T_numu)]

# Start the integration at a common temperature of 20 MeV, which corresponds to t ~ 2*10^{-3} s
T0 = 20 * 1e-3
t0 = 1./(2*Hubble(T0,T0,T0))

# Finish the calculation at t = 5*10^4 seconds, which will correspond to T ~ 5*10^{-3} MeV
t_max = 5e4

# Calculate
tvec = np.logspace(np.log10(t0),np.log10(t_max),num=300)
sol = odeint(dT_totdt, [T0,T0,T0], tvec, rtol = 1e-12, atol= 1e-12)

# Display Results
print "Neff        = ", round(Neff_func(sol[-1,0],sol[-1,1],sol[-1,2]),5)
print "Tgam/Tnu_e  = ", round(sol[-1,0]/sol[-1,1],5)
print "Tgam/Tnu_mu = ", round(sol[-1,0]/sol[-1,2],5)

# Calculate some Thermodynamic quantities
rho_vec, p_vec = np.zeros(len(tvec)),np.zeros(len(tvec))
for i in xrange(len(tvec)):
    rho_vec[i], p_vec[i]  = rho_tot(sol[i,0], sol[i,1], sol[i,2]), p_tot(sol[i,0], sol[i,1], sol[i,2])

# Store results
with open("SM/SM_emu.dat", 'w') as f:
    f.write("# SM with different Tnue and Tnumu = Tnutau "+"\n")
    f.write("# MB collision terms and QED corrections"+"\n")
    f.write("# Neff        = "),f.write("{:<12}".format(round(Neff_func(sol[-1,0],sol[-1,1],sol[-1,2]),5))),f.write("\n")
    f.write("# Tgam/Tnu_e  = "),f.write("{:<12}".format(round(sol[-1,0]/sol[-1,1],5))),f.write("\n")
    f.write("# Tgam/Tnu_mu = "),f.write("{:<12}".format(round(sol[-1,0]/sol[-1,2],5))),f.write("\n")

    f.write("{:<13}".format("#T_gamma (MeV)")), f.write("{:<13}".format("T_gam/T_nue")), f.write("{:<13}".format("T_gam/T_numu")),f.write("{:<14}".format("R_tot/T_gam^4 ")),f.write("{:<14}".format("P_tot/T_gam^4"+"\n"))
    for i in xrange(len(tvec)):
        f.write("{:.7E}".format(1e3*sol[i,0])),f.write("  "), f.write("{:<12}".format(round(sol[i,0]/sol[i,1],6))), f.write("{:<12}".format(round(sol[i,0]/sol[i,2],6))), f.write("{:.7E}".format(rho_vec[i]/sol[i,0]**4)),f.write("  "),f.write("{:.7E}".format(p_vec[i]/sol[i,0]**4)),f.write("\n")



