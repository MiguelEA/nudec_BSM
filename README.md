# NUDEC_BSM: Neutrino Decoupling Beyond the Standard Model

This code "NUDEC_BSM", has been developed by Miguel Escudero Abenza in order to solve for early Universe thermodynamics and neutrino decoupling following the simplified approach of ArXiv:1812.05605 [JCAP 1902 (2019) 007] and ArXiv:2001.XXXXX. If you use this code, please, cite these references. 

As of 10/01/2020:

There is a Mathematica and a Python version of NUDEC_BSM. The code consists of various scripts that calculate early Universe thermodynamics in various scenarios typically within the context of neutrino decoupling. 

The Python version contains the following scripts:

NUDEC_BSM.py          : This is a runner file that shows an example of how to run each of the models coded up in Python. 

nuDec_SM.py           : Solves for neutrino decoupling in the SM.

nuDec_SM_2_nu.py      : Solves for neutrino decoupling in the SM allowing evolving seperately  nu_e and nu_mu populations.

WIMP_e.py             : Solves for neutrino decoupling in the presence of a particle in thermal equilibrium with the electromagnetic component of the plasma.

WIMP_nu.py            : Solves for neutrino decoupling in the presence of a particle in thermal equilibrium with the neutrino component of the plasma.

WIMP_generic.py       : Solves for neutrino decoupling in the presence of a particle in thermal equilibrium with either the neutrino or electromagnetic component of the SM plasma, but which still interacts with the other component by means of annihilations. 

The header of each script contains the details on how to run and NUDEC_BSM.py contains an example for each case. 

The Mathematica version contains the following scripts:

BasicModules.nb     : Contains basic modules used in every model: QED finite temperature correctiosn, Thermodynamic formulae, SM interaction rates, constants, and parameters. When run it outputs BasicModules.m that can be loaded in any module.

Neff_SM.nb          : Solves for neutrino decoupling in the Standard Model. To run one should simply run the entire script and see the examples and output.

DarkRadiation.nb    : Solves for neutrino decoupling in the presence of Dark Radiation. One should simply run the entire script to find the thermodynamics. The only input parameter in this case is DNeff. 

NuScalar.nb         : Solves for the early Universe thermodynamics in the presence of a very light (eV<m<MeV) and weakly coupled (lambda < 10^{-9}) neutrinophilic scalar. There are two input parameters in this case Gamma_eff and mphi (MeV). In the particular scenario considered for this case the results have been shown to be very accurate for Gamma_eff > 10^{-3}. Note that for Gamma_eff < 10^{-3} the accuracy could be substantially lowered. 
