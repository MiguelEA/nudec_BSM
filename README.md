# nudec_BSM Neutrino Decoupling Beyond the Standard Model

UPDATE 14/08/19. 

Corrected the neutrino-neutrino and electron-neutrino reaction rates. Now the resulting Neff in the SM is 3.053. More details to appear in a new version of paper on the arxiv and in an erratum in JCAP.

Code has been updated to be compatible with both Python2 and Python3.

As of 13/12/18: 

This code "nudec_BSM", has been developed by Miguel Escudero Abenza in order to solve for the neutrino decoupling in the simplified approach to it developed in ArXiv:1812.05605, JCAP 1902 (2019) 007.

The user is referred to ArXiv:1812.05605 for the relevant equations and approximations that lead to them.

If you use this code, please cite: ArXiv:1812.05605.

The code "nudec_BSM" as of 13/12/18 contains various python scripts which evaluate N_eff as relevant for CMB observations, and provide the relevant thermodynamic variables for 20 MeV >~ T_\gamma >~ 0.005 MeV. They are:

nuDec_SM.py : Solves the neutrino decoupling in the SM

nuDec_SM_2_nu.py : Solves the neutrino decoupling in the SM allowing for a different temperature for the nu_e and nu_mu populations

WIMP_e.py : Solves the neutrino decoupling with a particle in thermal equilibrium with the electromagnetic component of the SM plasma.

WIMP_nu.py : Solves the neutrino decoupling with a particle in thermal equilibrium with the neutrino component of the SM plasma.

WIMP_generic.py : Solves the neutrino decoupling with a particle in thermal equilibrium with either the neutrino or electromagnetic component of the SM plasma, but which still interacts with the other component by means of annihilations.

The header of each script contains the details on how to run.
