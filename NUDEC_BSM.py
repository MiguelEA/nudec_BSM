# Last modified 27/05/2020. Miguel Escudero Abenza, miguel.escudero@kcl.ac.uk
# In this .py file we provide examples of how to solve neutrino decoupling in any of the models coded up in Python:
# 1) SM with a common neutrino temperature:             nuDec_SM.py
# 2) SM evolving separately nu_e and nu_mu-tau:         nuDec_SM_2nu.py
# 3) Thermal WIMP coupled to electrons:                 WIMP_e.py
# 4) Thermal WIMP coupled to neutrinos:                 WIMP_nu.py
# 5) Thermal WIMP annihilating differently to e and nu: WIMP_generic.py
import os
import numpy as np

# Calculate the SM Thermodynamics
print("\n"+"SM Thermodynamics:")
os.system("python nuDec_SM.py")

# Calculate the SM Thermodynamics evolving separately nu_e and nu_mu-tau
print("\n"+"SM Thermodynamics, Tnue != Tnumu :")
os.system("python nuDec_SM_2nu.py")

# Calculate the Thermodynamics with a thermal a Dark Photon of m = 5 MeV
MDM, gDM, SPIN = 5.0, 3.0, "B"
print("\n"+"Dark Photon m = 5 MeV:")
os.system("python WIMP_e.py "+"{:.2E}".format(MDM)+" "+str(gDM)+" "+SPIN+" aux_"+"{:.2E}".format(MDM)+"_"+str(gDM)+"_"+SPIN)

# Calculate the Thermodynamics with a thermal neutrinophilic Dirac fermion with m = 10 MeV
MDM, gDM, SPIN = 10.0, 4.0, "F"
print("\n"+"Neutrinophilic Dirac Fermion m = 10 MeV:")
os.system("python WIMP_nu.py "+"{:.2E}".format(MDM)+" "+str(gDM)+" "+SPIN+" aux_"+"{:.2E}".format(MDM)+"_"+str(gDM)+"_"+SPIN)

# Calculate the Thermodynamics with a thermal neutrinophilic Majoran fermion with m = 5 MeV
# with an annihilation cross section = 3*10^{-26} cm^3/s and annihilating in a ratio 4:1 to neutrinos:electrons
MDM, gDM, SPIN, type, BR = 7.0, 2.0, "F", "nu", 1./4.
print("\n"+"7 MeV Majorana WIMP with a branching ratio 4:1 to nu:e :")
os.system("python WIMP_generic.py "+"{:.2E}".format(MDM)+" "+str(gDM)+" "+SPIN+" "+type+" "+"{:.2E}".format(BR)+" aux_"+"{:.2E}".format(MDM)+"_"+str(gDM)+"_"+SPIN+"_"+"{:.2E}".format(BR))

