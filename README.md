# nudec_BSM The Neutrino Decoupling Beyond the Standard Model

This code "nudec_BSM", has been developed by Miguel Escudero Abenza in order to solve for the neutrino decoupling in the simplified approach to it developed in arxiv:1812.05605.

The user is referred to arxiv:1812.05605 for the relevant equations and approximations that lead to them.

If you use this code, please cite: arxiv:1812.05605

The code "nudec_BSM" as of 13/12/18 contains various python scripts which evaluate N_eff as relevant for CMB observations, and provide the relevant thermodynamic variables for 20 MeV >~ T_\gamma >~ 0.005 MeV. They are:

nuDec_SM.py : Solves the neutrino decoupling in the SM

nuDec_SM_2_nu.py : Solves the neutrino decoupling in the SM allowing for a different temperature for the nu_e and nu_mu populations

WIMP_e.py : Solves the neutrino decoupling with a particle in thermal equilibrium with the electromagnetic component of the SM plasma.

WIMP_nu.py : Solves the neutrino decoupling with a particle in thermal equilibrium with the neutrino component of the SM plasma.

WIMP_generic.py : Solves the neutrino decoupling with a particle in thermal equilibrium with either the neutrino or electromagnetic component of the SM plasma, but which still interacts with the other component by means of annihilations.

The header of each script contains the details on how to run.

The python version under which the code has been developed is 2.7. One can use Python 3 by making minor modifications in the print command and by changing xrange to range.
