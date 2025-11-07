# Neutrino Decoupling in and beyond the Standard Model

This code "nudec_BSM_v2" provides a fast and flexible description of neutrino decoupling in the early universe following the momentum-averaged approach. Independent Mathematica and Python implementations are provided, and should be easy to link with other codes, e.g. for BBN. The latest version has been developed by M. Escudero, G. Jackson, M. Laine and S. Sandner in [2511.XXXXX](https://arxiv.org/abs/2511.XXXXX) and builds on the earlier v1 from M. Escudero in [1812.05605](https://arxiv.org/abs/1812.05605) and [2001.04466](https://arxiv.org/abs/2001.04466). 

As of 07/11/2025:

New version 2 of the code being imported under branch ``nudec_v2``.

## Python implementation

This version is compatible with Python3 and contains the following scripts:

* ``test.py`` : Minimal example file that shows how to run the basic models in Python. 

* ``nudec.ipynb`` : Jupyter notebook walkthrough, including several plots. 

* ``BasicModules_source/nudec_v2.py`` : Main source file, includes classes needed to solve for neutrino decoupling. 


## Mathematica implementation

The Mathematica version contains the following scripts:

* ``Neff_SM.nb`` : Contains several illustrations of how to import and run the modules in Mathematica

* ``BasicModules_source/nudec_v2.m`` : Mathematica package file (see also notebook with same name).

## Pre-tabulated data

Stored under ``BasicModules_data/``. This includes electron mass dependent coefficients for neutrino interaction rates, as well as the interacting QED pressure and its first two derivatives. See sec.8 of [2511.XXXXX](https://arxiv.org/abs/2511.XXXXX).
