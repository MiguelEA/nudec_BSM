"""
  nudec_BSM_v2: neutrino decoupling in and beyond the Standard Model
  Copyright (C) 2025 M. Escudero, G. Jackson, M. Laine, S. Sandner

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version. See <https://www.gnu.org/licenses/>.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  [Usage:]

  example.py provides a minimal working example of how to evaluate the 
  temperature evolution as studied in the associated paper. A jupyter 
  notebook nudec.ipynb is also included The class is initialized as

    nudec   = nudec_source.NuDec([optional args below])

    ----------
    use_data_QED : [bool]
                   Indicates whether the QED equation of state is obtained 
                   from pre-tabulated files, e.g. like table III. The file 
                   directories are pointed to by the additional arguments, e.g. 

                     data_QED_p_int       -> QED_p_int.dat
                     data_QED_dp_dT_int   -> QED_dp_dT_int.dat
                     data_QED_d2p_dT2_int -> QED_d2p_dT2_int.dat

    Bessel       : [bool]
                   Switch for using expansions of the thermodynamic special 
                   functions in terms of K_n Bessel functions (appendix E). 
                   If set to False, these functions will be computed by 
                   numerical quadrature of their integral representations.

    interp       : [bool]
                   If True, will compute a grid of values for the thermodynamic 
                   functions and use that to construct an interpolations function 
                   on initialization of the class. This will allow for faster 
                   evaluation when solving the temperature evolution. The number 
                   of points on the grid is set by the additional argument:

                     interp_num     -> 1e4 (default)

    data_rates   : [string]
                   Points to the file which should be used for the energy and 
                   number density transfer rates, namely the default provided: 
                   -> rate_coefficients_neutrinos.dat (stored in data/)

  Once the class is ready, the following line will solve the necessary 
  initial value problem (from T_ini = 10 MeV to T_fin = 0.009 MeV by default). 
  The resulting solution (a vector of temperatures and chemical potentials) 
  can be accessed directly:

    [t, y]  = nudec.evolve([optional args below])

    ----------
    T_ini   : [float]
              Initial temperature in Mev. (default = 10 MeV) 

    T_fin   : [float]
              Final temperature in Mev. (default = 0.009 MeV) 

    mode    : [int]
              Indicates which 'mode' the evolution will run in (meaning which 
              system of ODEs will be solved). The options are as follows, 

                1 -> Tnu_e != Tnu_mu, mu_nue  = mu_numu = 0 [no oscillations & mu_nu  = 0]
                2 -> Tnu_e != Tnu_mu, mu_nue != mu_numu     [no oscillations & mu_nu != 0]
                3 -> Tnu_e  = Tnu_mu, mu_nue  = mu_numu = 0 [w/ oscillations & mu_nu  = 0]
                4 -> Tnu_e  = Tnu_mu, mu_nue  = mu_numu     [w/ oscillations & mu_nu != 0]

    order   : [int]
              Selection of the perturbative order for the (interacting) QED 
              equation of state. Options are as follows, 

                0, 1 -> O(e^0)
                2    -> O(e^2,non-log term)
                3    -> O(e^3)
                4    -> O(e^4)
                5    -> O(e^5)
           {and 6    -> O(e^5+e^2 log term), but only if use_data_QED = True!}

    me      : [float]
              Electron mass in MeV. (default stored in NuDec_Const) 

    mpl     : [float]
              (reduced) Planck mass in MeV. (default stored in NuDec_Const) 

    GF      : [float]
              Fermi coupling constant in MeV^{-2}. (default stored in NuDec_Const) 

    e       : [float]
              Electromagnetic coupling constant e=sqrt{4.pi.alpha}. (default stored in NuDec_Const) 

    geL, geR, gmuL, gmuR : [floats]
              L/R lepton couplings.

    info    : [bool]
              If set to True, will print out intermediate information.

    ode_rtol, ode_atol: [floats]
              Sets the relative and absolute tolerances of the python IVP solver, 
              for method='LSODA'. (default stored in NuDec_Const)


  Physical constants are stored in the 'NuDec_Const' class, along with 
  several other parameters to control the numerical tolerances, 
  quadrature limits, and series expansions. 

  [Notes:]

  Some comments below refer to specific places in ArXiv:2511.XXXXX

"""
import os
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.integrate import quad, solve_ivp
from scipy.special import kn, zeta, spence, kv
from scipy.integrate import quad
import warnings

class NuDec_Const:
    me      = 0.510999              # electron mass in MeV the SM
    GF      = 1.1663787e-11         # Fermi constant in MeV^{-2}
    mpl     = 1.220890e22           # Planck mass in MeV
    alpha   = 1./137.035999084      # fine-structure constant (alpha=e^2/4.pi)
    e       = np.sqrt(4.*np.pi*alpha)

    m_nu        = 0     # assume massless neutrinos
    g_nu        = 2     # 2 internal degrees of freedom for neutrino + anti-neutrino

    T_ini       = 10.0
    T_fin       = 0.009

    geL, geR, gmuL, gmuR = 0.727, 0.233, -0.273, 0.233 # left and right nu-e couplings as relevant for E < 10 MeV

    MeVtoSec    = 1/(6.58212e-22)       # conversion factor to transform MeV^-1 into seconds

    # Temperature of the CMB today
    T_cmb = 2.7255*(1./11604.51812) # eV
    # critical energy density today
    rho_c = 8.0959*1e-11            # eV^4

    # QED perturbation order: 0, 1 -> O(e^0) ; 2 -> O(e^2,non-log term) ; 3 -> O(e^3) ; 4 -> O(e^4) ; 5 -> O(e^5)
    # [and 6 -> O(e^5+e^2 log term), but only if use_data_QED = True!]
    order       = 3

    # option to evolve with or without chemical potential and e \neq mu,tau temp
    mode        = 4

    # truncation order of the Bessel expansion for the special functions...
    # ...used to calculate thermodynamical quantities
    Bessel_Max  = 30

    # settings for the ODE solve accuracy and precision
    ode_rtol    = 1e-8
    ode_atol    = 1e-8

    # numerical integration accuracy and precision used for special function in the thermodynamics class
    quad_rtol   = 1e-6
    quad_atol   = 1e-6

    # should intermediate results/information be printed?
    info        = True

    # avoid numerical problems for tau -> 0 in the QED pressure
    tau_min_QED_Pressure    = 1e-5

    # tau cut-off for special functions
    tau_min_sf_series = np.max([0.2,10/Bessel_Max])
    tau_max_sf_series = 30

    # set path to current working directory 
    cwd         = os.getcwd()




class NuDec:
    def __init__(self, data_rates: str= NuDec_Const.cwd + '/data/rate_coefficients_neutrinos.dat', 
                    use_data_QED: bool= False, 
                    data_QED_p_int: str= NuDec_Const.cwd + '/data/QED_p_int.dat',
                    data_QED_dp_dT_int: str= NuDec_Const.cwd + '/data/QED_dp_dT_int.dat',
                    data_QED_d2p_dT2_int: str= NuDec_Const.cwd + '/data/QED_d2p_dT2_int.dat', 
                    Bessel: bool= True, interp: bool= True, interp_num: int= 1e4) -> None:
        self.Bessel         = Bessel
        self.use_data_QED   = use_data_QED
        if Bessel and use_data_QED:
            warnings.warn("In class NuDec: Arguments Bessel and use_data_QED cannot be simultaneously true.\n" \
                            "-> Code continues assuming Bessel= False, use_data_QED= True.", category=UserWarning, stacklevel=2)
        if use_data_QED:
            print("Initializing NuDec class with: use_data_QED= True.")
            print("Electron mass and QED coupling are fixed to their SM values and cannot be changed.")

        self.interp         = interp
        self.thermo         = self.Thermo(use_data_QED, data_QED_p_int, data_QED_dp_dT_int, data_QED_d2p_dT2_int, 
                                          Bessel, interp, interp_num)
        self.rates          = self.Rates(data_rates) 



    # Evolution equation for the case: Tnu_e != Tnu_mu & mu_nue = mu_numu = 0
    # "no oscillations and mu_nu = 0"
    def ode_fun_mode1(self, y: list, me: float= NuDec_Const.me, order: int= NuDec_Const.order, 
                      mpl: float= NuDec_Const.mpl, GF: float= NuDec_Const.GF, 
                      geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, 
                      gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, e: float= NuDec_Const.e) -> float:
        T_gam, T_nue, T_numu, z  = y #unpack
        mu_nu       = 0

        H           = self.Hubble(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, order=order, mpl=mpl, e=e)

        Rho_QED    = self.thermo.Rho_QED_int(T=T_gam, m=me, order=order, e=e)
        P_QED      = self.thermo.P_QED_int(T=T_gam, m=me, order=order, e=e)
        dRhodt_QED  = 2.*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF) \
                    +self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, geL=geL, geR=geR, GF=GF)

        Rho_e       = self.thermo.Rho_FD(T=T_nue, mu=mu_nu)
        dRhodT_e    = 3.*self.thermo.dP_dT_nu(T=T_nue, mu=mu_nu)
        dRhodt_e    = self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, geL=geL, geR=geR, GF=GF)

        Rho_mu       = 2*self.thermo.Rho_FD(T=T_numu, mu=mu_nu)
        dRhodT_mu    = 2*3.*self.thermo.dP_dT_nu(T=T_numu, mu=mu_nu)
        dRhodt_mu    = 2*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF)


        dTgam_dt    = (-3.*H*( Rho_QED + P_QED) - dRhodt_QED) / (T_gam*self.thermo.d2P_dT2_QED(T=T_gam, m=me, order=order, e=e)) 
        dTnue_dt    = (-4.*H*Rho_e + dRhodt_e)/(dRhodT_e)
        dTnumu_dt   = (-4.*H*Rho_mu + dRhodt_mu)/(dRhodT_mu)
        dz_dt       = z*(H + dTgam_dt/T_gam)


        return [dTgam_dt, dTnue_dt, dTnumu_dt, dz_dt] 

    # Evolution equation for the case: Tnu_e != Tnu_mu & mu_nue != mu_numu
    # "no oscillations and mu_nu != 0"
    def ode_fun_mode2(self, y: list, me: float= NuDec_Const.me, order: int= NuDec_Const.order, 
                      mpl: float= NuDec_Const.mpl, GF: float= NuDec_Const.GF, 
                      geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, 
                      gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, e: float= NuDec_Const.e) -> float:

        T_gam, T_nue, T_numu, mu_nue, mu_numu, z  = y #unpack

        H           = self.Hubble(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, order=order, mpl=mpl, e=e)

        Rho_QED     = self.thermo.Rho_QED_int(T=T_gam, m=me, order=order, e=e)
        P_QED       = self.thermo.P_QED_int(T=T_gam, m=me, order=order, e=e)
        dRhodt_QED  = 2.*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF)\
                    +self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, geL=geL, geR=geR, GF=GF)


        n_e         = self.thermo.n_FD(T=T_nue, mu=mu_nue, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu)
        dndmu_e     = self.thermo.dn_dmu_nu(T=T_nue, mu=mu_nue)
        dndT_e      = self.thermo.dn_dT_nu(T=T_nue, mu=mu_nue)
        Rho_e       = self.thermo.Rho_FD(T=T_nue, mu=mu_nue)
        P_e         = Rho_e/3.
        dRhodmu_e   = 3.*self.thermo.dP_dmu_nu(T=T_nue, mu=mu_nue)
        dRhodT_e    = 3.*self.thermo.dP_dT_nu(T=T_nue, mu=mu_nue)
        dRhodt_e    = self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, geL=geL, geR=geR, GF=GF)
        dndt_e      = self.rates.DeltaN_nue(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, geL=geL, geR=geR, GF=GF)

        n_mu        = 2.*self.thermo.n_FD(T=T_numu, mu=mu_numu, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu)
        dndmu_mu    = 2.*self.thermo.dn_dmu_nu(T=T_numu, mu=mu_numu)
        dndT_mu     = 2.*self.thermo.dn_dT_nu(T=T_numu, mu=mu_numu)
        Rho_mu      = 2.*self.thermo.Rho_FD(T=T_numu, mu=mu_numu)
        P_mu        = Rho_mu/3.
        dRhodmu_mu  = 2.*3.*self.thermo.dP_dmu_nu(T=T_numu, mu=mu_numu)
        dRhodT_mu   = 2.*3.*self.thermo.dP_dT_nu(T=T_numu, mu=mu_numu)
        dRhodt_mu   = 2.*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF)
        dndt_mu     = 2.*self.rates.DeltaN_numu(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF)



        dTgam_dt    = (-3.*H*( Rho_QED + P_QED) - dRhodt_QED) / (T_gam*self.thermo.d2P_dT2_QED(T=T_gam, m=me, order=order, e=e))

        dTnue_dt    = (-3.*H*( (Rho_e+P_e)*dndmu_e - n_e*dRhodmu_e) + dndmu_e*dRhodt_e - dRhodmu_e*dndt_e ) /(dndmu_e*dRhodT_e - dndT_e*dRhodmu_e)
        dTnumu_dt   = (-3.*H*( (Rho_mu+P_mu)*dndmu_mu - n_mu*dRhodmu_mu) + dndmu_mu*dRhodt_mu - dRhodmu_mu*dndt_mu )/(dndmu_mu*dRhodT_mu - dndT_mu*dRhodmu_mu)

        dmue_dt     = -(-3.*H*( (Rho_e+P_e)*dndT_e - n_e*dRhodT_e) + dndT_e*dRhodt_e - dRhodT_e*dndt_e )/(dndmu_e*dRhodT_e - dndT_e*dRhodmu_e)
        dmumu_dt    = -(-3.*H*( (Rho_mu+P_mu)*dndT_mu - n_mu*dRhodT_mu) + dndT_mu*dRhodt_mu - dRhodT_mu*dndt_mu )/(dndmu_mu*dRhodT_mu - dndT_mu*dRhodmu_mu)

        dz_dt       = z*(H + dTgam_dt/T_gam)

        return [dTgam_dt, dTnue_dt, dTnumu_dt, dmue_dt, dmumu_dt, dz_dt] 

    # Evolution equation for the case: Tnu_e = Tnu_mu & mu_nue = mu_numu = 0
    # "with oscillations and mu_nu = 0"
    def ode_fun_mode3(self, y: list, me: float= NuDec_Const.me, order: int= NuDec_Const.order, 
                      mpl: float= NuDec_Const.mpl, GF: float= NuDec_Const.GF, 
                      geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, 
                      gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, e: float= NuDec_Const.e) -> float:

        T_gam, T_nu, z  = y  # unpack array
        mu_nu       = 0

        H           = self.Hubble(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, order=order, mpl=mpl, e=e)

        Rho_QED    = self.thermo.Rho_QED_int(T=T_gam, m=me, order=order, e=e)
        P_QED      = self.thermo.P_QED_int(T=T_gam, m=me, order=order, e=e)


        Rho         = 3.*self.thermo.Rho_FD(T=T_nu, mu=mu_nu)
        P           = Rho/3.
        dRhodT      = 9.*self.thermo.dP_dT_nu(T=T_nu, mu=mu_nu)

        dRhodt      = 2.*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF) \
                    +self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, geL=geL, geR=geR, GF=GF)


        dTgam_dt    = (-3.*H*( Rho_QED + P_QED) - dRhodt) / (T_gam*self.thermo.d2P_dT2_QED(T=T_gam, m=me, order=order, e=e)) 
        dTnu_dt     = - H*T_nu + dRhodt/dRhodT
        dz_dt       = z*(H + dTgam_dt/T_gam)



        return [dTgam_dt, dTnu_dt, dz_dt] 

    # Evolution equation for the case: Tnu_e = Tnu_mu & mu_nue = mu_numu
    # "with oscillations and mu_nu != 0"
    def ode_fun_mode4(self, y: list, me: float= NuDec_Const.me, order: int= NuDec_Const.order, 
                      mpl: float= NuDec_Const.mpl, GF: float= NuDec_Const.GF, 
                      geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, 
                      gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, e: float= NuDec_Const.e) -> float:

        T_gam, T_nu, mu_nu, z  = y  # unpack array

        H           = self.Hubble(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, order=order, mpl=mpl, e=e)

        Rho_QED    = self.thermo.Rho_QED_int(T=T_gam, m=me, order=order, e=e)
        P_QED      = self.thermo.P_QED_int(T=T_gam, m=me, order=order, e=e)


        n           = 3.*self.thermo.n_FD(T=T_nu, mu=mu_nu, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu)
        dndmu       = 3.*self.thermo.dn_dmu_nu(T=T_nu, mu=mu_nu)
        dndT        = 3.*self.thermo.dn_dT_nu(T=T_nu, mu=mu_nu)
        Rho         = 3.*self.thermo.Rho_FD(T=T_nu, mu=mu_nu)
        P           = Rho/3.
        dRhodmu     = 9.*self.thermo.dP_dmu_nu(T=T_nu, mu=mu_nu)
        dRhodT      = 9.*self.thermo.dP_dT_nu(T=T_nu, mu=mu_nu)

        dRhodt      = 2.*self.rates.DeltaRho_numu(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF) \
                    +self.rates.DeltaRho_nue(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, geL=geL, geR=geR, GF=GF)
        dndt        = 2.*self.rates.DeltaN_numu(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, gmuL=gmuL, gmuR=gmuR, GF=GF) \
                    +self.rates.DeltaN_nue(T_gam=T_gam, T_nue=T_nu, T_numu=T_nu, mu_nue=mu_nu, mu_numu=mu_nu, me=me, geL=geL, geR=geR, GF=GF)



        dTgam_dt    = (-3.*H*( Rho_QED + P_QED ) - dRhodt) / (T_gam*self.thermo.d2P_dT2_QED(T=T_gam, m=me, order=order, e=e))
        dTnu_dt     = (-3.*H*( (Rho+P)*dndmu - n*dRhodmu )  + dndmu*dRhodt - dRhodmu*dndt ) / (dndmu*dRhodT - dndT*dRhodmu )
        dmu_dt      = - (-3.*H*( (Rho+P)*dndT - n*dRhodT )  + dndT*dRhodt - dRhodT*dndt) / (dndmu*dRhodT - dndT*dRhodmu )
        dz_dt       = z*(H + dTgam_dt/T_gam)

        return [dTgam_dt, dTnu_dt, dmu_dt, dz_dt] 



    def evolve(self, T_ini: float= NuDec_Const.T_ini, T_fin: float= NuDec_Const.T_fin, me: float= NuDec_Const.me, 
                order: int= NuDec_Const.order, mpl: float= NuDec_Const.mpl,
                GF: float= NuDec_Const.GF, geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, 
                gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, e: float= NuDec_Const.e,
                mode: int= NuDec_Const.mode, info: bool= NuDec_Const.info, 
                ode_rtol: float= NuDec_Const.ode_rtol, ode_atol: float= NuDec_Const.ode_atol):

        if mode not in (1, 2, 3, 4):
            raise ValueError("Argument 'mode' must be 1, 2, 3 or 4.")
        if order not in (0, 1, 2, 3, 4, 5, 6):
            warnings.warn("Argument 'order' must be either 0, 1, 2, 3, 4, 5, 6 " \
                            "Code re-sets to default order = 3.", category=UserWarning, stacklevel=2)
            order = NuDec_Const.order


        t_ini = 1./(2*self.Hubble(T_gam=T_ini, T_nue=T_ini, T_numu=T_ini, mu_nue=0, mu_numu=0, me=me, order=order, mpl=mpl, e=e))
        t_fin = 1./(2*self.Hubble(T_gam=T_fin, T_nue=T_fin/1.4, T_numu=T_fin/1.4, mu_nue=0, mu_numu=0, me=me, order=order, mpl=mpl, e=e)) # estimated

        if info:
            print(f"\n solve system from T_ini = {T_ini:.3f} to T_fin = {T_fin:.3f} [MeV]")
            print(f" --> t_ini = {t_ini:.2E} [s] to t_fin = {t_fin:.2E} [s]\n")

        if mode==1:
            y_ini   = [T_ini,T_ini,T_ini,1]
        elif mode==2:
            y_ini   = [T_ini,T_ini,T_ini,0,0,1]
        elif mode==3:
            y_ini   = [T_ini,T_ini,1]
        else:
            y_ini   = [T_ini,T_ini,0,1] 


        def dT_totdt(t, y):
            if info:
                progress = 100 * (t - t_ini) / (t_fin - t_ini)
                print(f"\rProgress: {progress:.3f} %", end= "", flush= True)
            if mode==1:
                return self.ode_fun_mode1(y=y, me=me, order=order, mpl=mpl, GF=GF, geL=geL, geR=geR, gmuL=gmuL, gmuR=gmuR, e=e)
            elif mode==2:
                return self.ode_fun_mode2(y=y, me=me, order=order, mpl=mpl, GF=GF, geL=geL, geR=geR, gmuL=gmuL, gmuR=gmuR, e=e)
            elif mode==3:
                return self.ode_fun_mode3(y=y, me=me, order=order, mpl=mpl, GF=GF, geL=geL, geR=geR, gmuL=gmuL, gmuR=gmuR, e=e)
            else:
                return self.ode_fun_mode4(y=y, me=me, order=order, mpl=mpl, GF=GF, geL=geL, geR=geR, gmuL=gmuL, gmuR=gmuR, e=e)


        sol     = solve_ivp( dT_totdt, t_span=[t_ini, t_fin], y0=y_ini,
                method='LSODA', rtol= ode_rtol, atol= ode_atol)
        t, y    = np.array(sol.t), np.array(sol.y)


        # Output formating...
        if mode==1:
            idx_munu    = 2
            zero_row    = np.zeros((1, y.shape[1]), dtype=y.dtype)
            y           = np.vstack((y[:idx_munu+1], zero_row, zero_row, y[idx_munu+1:]))
        if mode==3:
            # for output add 2 munu rows and 1 copied Tnu row
            idx_Tnu     = 1
            Tnu_copy    = y[[idx_Tnu], :]
            zero_row    = np.zeros((1, y.shape[1]), dtype=y.dtype)
            y           = np.vstack((y[:idx_Tnu+1], Tnu_copy, zero_row, zero_row, y[idx_Tnu+1:]))   
        if mode==4:
            idx_Tnu, idx_munu     = 1, 2
            Tnu_copy, munu_copy = y[[idx_Tnu], :], y[[idx_munu], :] 
            y   = np.vstack((y[:idx_Tnu+1, :], Tnu_copy, y[idx_Tnu+1:idx_munu+1, :], munu_copy, y[idx_munu+1:, :]))


        return [t, y] # = [t, [T_gam, T_nue, T_numu, mu_nue, mu_numu, z]]



    #   ------------------------
    #   observables of interest
    #   ------------------------

    def Neff(self, T_gam, T_nue, T_numu, mu_nue, mu_numu, me: float= NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
        """
        effective number of neutrino species: cf. eq.(5.2)
        """
        prefac  =  8./7.*(11./4.)**(4./3.)

        return prefac*(self.thermo.Rho_FD(T= T_nue, mu= mu_nue, m= NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu) + 2*self.thermo.Rho_FD(T= T_numu, mu= mu_numu, m= NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu))/self.thermo.Rho_QED_int(T= T_gam, m= me, order= order, e= e)

    def gstar_rho(self, T_gam, T_nue, T_numu, mu_nue, mu_numu, me: float= NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
        """
        effective d.o.f. from energy density: cf. eq.(5.1)
        """
        prefac = 30./np.pi**2

        return prefac*( self.thermo.Rho_tot(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, order=order, e=e) )/T_gam**4

    def gstar_s(self, T_gam, T_nue, T_numu, mu_nue, mu_numu, me: float= NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
        """
        kinetic theory definition of entropy: cf. eq.(C.3) 
        """
        prefac  =  45./2./np.pi**2 

        return prefac*(self.thermo.s_FD(T=T_nue, mu=mu_nue, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu) + 2*self.thermo.s_FD(T=T_numu, mu=mu_numu, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu) + self.thermo.s_QED_int(T=T_gam, m=me, order=order, e=e))/T_gam**3

    def m_over_Omega(self, T_gam, T_nue, T_numu, mu_nue, mu_numu) -> float:
        """
        non-relativistic energy density: cf. eq.(C.11)
        """

        return 1./( ( self.thermo.n_FD(T= T_nue, mu= mu_nue, m= NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu) + 2*self.thermo.n_FD(T= T_numu, mu= mu_numu, m= NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu) )*NuDec_Const.T_cmb**3/(3*T_gam**3*NuDec_Const.rho_c) )

    def heff(self, z: float):
        """
        non-equilibrium entropy coefficient: cf. eq.(5.6) and eq.(D.15)
        """
        T       = NuDec_Const.T_ini
        s_ini   = self.thermo.s_QED_int(T=T, m=NuDec_Const.me, order=NuDec_Const.order, e=NuDec_Const.e)
        h0ini   = (45./np.pi**2)*( .5*s_ini/T**3 ) + 21./4.

        return h0ini/z**3


    # Hubble expansion rate:
    def Hubble(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, order: int= NuDec_Const.order, mpl: float= NuDec_Const.mpl, e: float= NuDec_Const.e) -> float:
        """ 
        expansion rate from energy density: cf. eq.(2.7)
        """
        return NuDec_Const.MeVtoSec*np.sqrt(self.thermo.Rho_tot(T_gam=T_gam, T_nue=T_nue, T_numu=T_numu, mu_nue=mu_nue, mu_numu=mu_numu, me=me, order=order, e=e)*8*np.pi/(3*mpl**2))




    class Rates:

        def __init__(self, data_rates) -> None:
            # improved MB rate coefficients (cf. sec. 2.B)
            # read in tabulation of mass correction factors: [ tau = 0,...,50 ]
            self.data       = np.loadtxt(data_rates,comments='#')
            tau             = self.data[:, 0]
            ycols           = self.data[:, range(1,13)].T  

            [self.f_a1_int, self.f_a2_int, self.f_a3_int, self.f_a4_int] = [make_interp_spline(tau, ycols[i]) for i in range(4)]
            [self.f_s1_int, self.f_s2_int, self.f_s3_int, self.f_s4_int] = [make_interp_spline(tau, ycols[i+4]) for i in range(4)]
            [self.f_n1_int, self.f_n2_int, self.f_n3_int, self.f_n4_int] = [make_interp_spline(tau, ycols[i+8]) for i in range(4)]

            # zero mass limits:
            [self.f_a1_0, self.f_s1_0, self.f_n1_0, self.f_a3_0, self.f_n3_0] = [ycols[i][0] for i in [0,4,8,2,10]]




        # energy & number density transfer rates: cf. eqs.(2.22)-(2.26)
        # omit _F_s_34, because f_{s3} and f_{s4} vanish

        def F_a_12(self, T1: float, T2: float, mu1: float = 0, mu2: float= 0) -> float:
            return T1**9*np.exp(2.*mu1/T1) - T2**9*np.exp(2.*mu2/T2)

        def F_a_34(self, T1: float, T2: float, mu1: float = 0, mu2: float= 0) -> float:
            return ( (T1*T2)**(4.5) )*( np.exp(2.*mu2/T2) - np.exp(2.*mu1/T1) )

        def F_s_12(self, T1: float, T2: float, mu1: float = 0, mu2: float= 0) -> float:
            return T1**4*T2**4*(T1-T2)*np.exp(mu1/T1)*np.exp(mu2/T2)

        def F_n_12(self, T1: float, T2: float, mu1: float = 0, mu2: float= 0) -> float:
            return ( T1**8*np.exp(2.*mu1/T1) - T2**8*np.exp(2.*mu2/T2) )

        def F_n_34(self, T1: float, T2: float, mu1: float = 0, mu2: float= 0) -> float:
            return (T1*T2)**4*( np.exp(2.*mu2/T2) - np.exp(2.*mu1/T1) )



        def DeltaRho_nue(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, GF: float= NuDec_Const.GF) -> float:

            """
            energy exchange into the nu_e ensemble: cf. eqs.(2.27)-(2.28) [alpha->e]
            """
            X = 4*(geL**2 + geR**2)
            Y = 4*(geL*geR) # ~ V^2-A^2 coupling

            # first, nu-nubar production/annihilation from QED plasma:
            fa1 = self.f_a1_int(me/T_gam)
            fa2 = self.f_a2_int(me/T_gam)
            fa3 = self.f_a3_int(me/T_gam)
            fa4 = self.f_a4_int(me/T_gam)
            res1  = ( X*fa1 + Y*fa2 )*self.F_a_12(T_gam,T_nue,0,mu_nue)\
                    + ( X*fa3 + Y*fa4 )*self.F_a_34(T_gam,T_nue,0,mu_nue)

            # then, nu & nubar scattering on the QED plasma:
            fs1 = self.f_s1_int(me/T_gam)
            fs2 = self.f_s2_int(me/T_gam)
            res2  = 7./4.*( X*fs1 + Y*fs2 )*self.F_s_12(T_gam,T_nue,0,mu_nue)

            # then, same exchanges with the nu_mu ensemble: (factor 2 for flavour)
            res3    = 2*( self.f_a1_0*self.F_a_12(T_numu,T_nue,mu_numu,mu_nue)\
                        + self.f_a3_0*self.F_a_34(T_nue,T_nue,mu_numu,mu_nue) \
                        + (7./4.)*self.f_s1_0*self.F_s_12(T_numu,T_nue,mu_numu,mu_nue) )

            return 32.*NuDec_Const.MeVtoSec*GF**2/np.pi**5 * (res1 + res2 + res3)



        def DeltaRho_numu(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, GF: float= NuDec_Const.GF) -> float:
            """
            energy exchange into the nu_mu ensemble: cf. eqs.(2.27)-(2.28) [alpha->mu]
            """
            X = 4*(gmuL**2 + gmuR**2)
            Y = 4*(gmuL*gmuR) # ~ V^2-A^2 coupling

            # first, nu-nubar production/annihilation from QED plasma:
            fa1 = self.f_a1_int(me/T_gam)
            fa2 = self.f_a2_int(me/T_gam)
            fa3 = self.f_a3_int(me/T_gam)
            fa4 = self.f_a4_int(me/T_gam)
            res1  = ( X*fa1 + Y*fa2 )*self.F_a_12(T_gam,T_numu,0,mu_numu)\
                  + ( X*fa3 + Y*fa4 )*self.F_a_34(T_gam,T_numu,0,mu_numu)

            # then, nu & nubar scattering on the QED plasma:
            fs1 = self.f_s1_int(me/T_gam)
            fs2 = self.f_s2_int(me/T_gam)
            res2  = 7./4. * ( X*fs1 + Y*fs2 )*self.F_s_12(T_gam,T_numu,0,mu_numu)

            # then, same exchanges with the nu_mu ensemble: (note the order of arguments)
            res3    = - ( self.f_a1_0*self.F_a_12(T_numu,T_nue,mu_numu,mu_nue)\
                        + self.f_a3_0*self.F_a_34(T_nue,T_nue,mu_numu,mu_nue) \
                        + (7./4.)*self.f_s1_0*self.F_s_12(T_numu,T_nue,mu_numu,mu_nue) )

            return 32.*NuDec_Const.MeVtoSec*GF**2/np.pi**5 * (res1 + res2 + res3)

        def DeltaN_nue(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, geL: float= NuDec_Const.geL, geR: float= NuDec_Const.geR, GF: float= NuDec_Const.GF) -> float:

            """
            number density exchange into the nu_e ensemble: cf. eqs.(2.29)-(2.30) [alpha->e]
            """
            X = 4*(geL**2 + geR**2)
            Y = 4*(geL*geR) # ~ V^2-A^2 coupling

            # first, nu-nubar production/annihilation from QED plasma:
            fn1 = self.f_n1_int(me/T_gam)
            fn2 = self.f_n2_int(me/T_gam)
            fn3 = self.f_n3_int(me/T_gam)
            fn4 = self.f_n4_int(me/T_gam)
            res1  = ( X*fn1 + Y*fn2 )*self.F_n_12(T_gam,T_nue,0,mu_nue)\
                + ( X*fn3 + Y*fn4 )*self.F_n_34(T_gam,T_nue,0,mu_nue)

            # then, same exchanges with the nu_mu ensemble: (factor 2 for flavour)
            res2    = 2*( self.f_n1_0*self.F_n_12(T_numu,T_nue,mu_numu,mu_nue)\
                        + self.f_n3_0*self.F_n_34(T_nue,T_nue,mu_numu,mu_nue))

            return 8.*NuDec_Const.MeVtoSec*GF**2/np.pi**5 * ( res1 + res2 )



        def DeltaN_numu(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, gmuL: float= NuDec_Const.gmuL, gmuR: float= NuDec_Const.gmuR, GF: float= NuDec_Const.GF) -> float:
            """
            number density exchange into the nu_mu ensemble: cf. eqs.(2.29)-(2.30) [alpha->mu]
            """ 

            X = 4*(gmuL**2 + gmuR**2)
            Y = 4*(gmuL*gmuR) # ~ V^2-A^2 coupling

            # first, nu-nubar production/annihilation from QED plasma:
            fn1 = self.f_n1_int(me/T_gam)
            fn2 = self.f_n2_int(me/T_gam)
            fn3 = self.f_n3_int(me/T_gam)
            fn4 = self.f_n4_int(me/T_gam)
            res1  = ( X*fn1 + Y*fn2 )*self.F_n_12(T_gam,T_numu,0,mu_numu)\
                  + ( X*fn3 + Y*fn4 )*self.F_n_34(T_gam,T_numu,0,mu_numu)

            # then, same exchanges with the nu_mu ensemble: (note the order of arguments)
            res2    = - ( self.f_n1_0*self.F_n_12(T_numu,T_nue,mu_numu,mu_nue)\
                        + self.f_n3_0*self.F_n_34(T_nue,T_nue,mu_numu,mu_nue))

            return 8.*NuDec_Const.MeVtoSec*GF**2/np.pi**5 * ( res1 + res2 )




    class Thermo:

        def __init__(self, use_data_QED,
                    data_QED_p_int,
                    data_QED_dp_dT_int,
                    data_QED_d2p_dT2_int,
                    Bessel, interp, interp_num) -> None:

            self.c4, self.c5    = -0.0015941171507113735422, -0.0006971656122053014132 # eq.(D.41), eq.(D.42)

            self.Bessel         = Bessel
            self.use_data_QED   = use_data_QED
            if Bessel and use_data_QED:
                warnings.warn("In class Thermo: Arguments Bessel and use_data_QED cannot be simultaneously true.\n" \
                                "-> Code continues assuming Bessel= False, use_data_QED= True.", category=UserWarning, stacklevel=2)

            self.interp     = interp
            if use_data_QED==False:
                self.interp_num = interp_num
                self.sf         = self.SpecialFunctions(Bessel= Bessel, interp= interp, interp_num= interp_num)
            else:
                self.p_dat      = np.loadtxt(data_QED_p_int, comments='#')
                self.dp_dat     = np.loadtxt(data_QED_dp_dT_int, comments='#')
                self.d2p_dat    = np.loadtxt(data_QED_d2p_dT2_int, comments='#')

                self.T_gam_dat  = self.p_dat[:,0]
                self.p_cols     = self.p_dat[:,range(1,7)].T
                self.dp_cols    = self.dp_dat[:,range(1,7)].T
                self.d2p_cols   = self.d2p_dat[:,range(1,7)].T

                print(f"\rInterpolating QED data file ...", end= "", flush= True)
                self.p_int_0, self.p_int_2_non_ln, self.p_int_2_ln, self.p_int_3, self.p_int_4, self.p_int_5 = [make_interp_spline(self.T_gam_dat, self.p_cols[i]) for i in range(6)]
                self.dp_int_0, self.dp_int_2_non_ln, self.dp_int_2_ln, self.dp_int_3, self.dp_int_4, self.dp_int_5 = [make_interp_spline(self.T_gam_dat, self.dp_cols[i]) for i in range(6)]
                self.d2p_int_0, self.d2p_int_2_non_ln, self.d2p_int_2_ln, self.d2p_int_3, self.d2p_int_4, self.d2p_int_5 = [make_interp_spline(self.T_gam_dat, self.d2p_cols[i]) for i in range(6)]
                print(f"\rInterpolating QED data file ... done", end= "", flush= True)

            self.polylog        = self.PolyLog()


        class SpecialFunctions:

            def __init__(self, Bessel, interp, interp_num) -> None:
                self.Bessel = Bessel

                if interp:
                    self.array_tau  = np.logspace(-6, 6, num= int(interp_num))

                    print(f"\rj interpolation ...                                           ", end= "", flush= True)
                    self.array_j        = [self.j(i, False) for i in self.array_tau]
                    self.j_interp       = make_interp_spline(self.array_tau, self.array_j)
                    print(f"\rj interpolation ... done  ", end= "", flush= True)

                    print(f"\rjp interpolation ...                                           ", end= "", flush= True)
                    self.array_jp        = [self.jp(i, False) for i in self.array_tau]
                    self.jp_interp       = make_interp_spline(self.array_tau, self.array_jp)
                    print(f"\rjp interpolation ... done  ", end= "", flush= True)

                    print(f"\rJ interpolation ...                                           ", end= "", flush= True)
                    self.array_J        = [self.J(i, False) for i in self.array_tau]
                    self.J_interp       = make_interp_spline(self.array_tau, self.array_J)
                    print(f"\rJ interpolation ... done  ", end= "", flush= True)

                    print(f"\rY interpolation ...                                           ", end= "", flush= True)
                    self.array_Y        = [self.Y(i, False) for i in self.array_tau]
                    self.Y_interp       = make_interp_spline(self.array_tau, self.array_Y)
                    print(f"\rY interpolation ... done  ", end= "", flush= True)

                    print(f"\rk interpolation ...                                           ", end= "", flush= True)
                    self.array_k        = [self.k(i, False) for i in self.array_tau]
                    self.k_interp       = make_interp_spline(self.array_tau, self.array_k)
                    print(f"\rk interpolation ... done  ", end= "", flush= True)

                    print(f"\rK interpolation ...                                           ", end= "", flush= True)
                    self.array_K        = [self.K(i, False) for i in self.array_tau]
                    self.K_interp       = make_interp_spline(self.array_tau, self.array_K)
                    print(f"\rK interpolation ... done  ", end= "", flush= True)


                    print(f"\rFinished constructing interpolation functions.        ", end= "", flush= True)


            def j(self, tau: float, interp: bool= True) -> float:
                """
                cf. eq.(D.16) and eq.(E.12)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series: 
                    zetapminus2 =-0.0304484570583933
                    return ( .5 + tau**2*(7/2.)*(zetapminus2) )/np.pi**2
                else:
                    if interp:
                        return self.j_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1)**i*(1+i)*tau*kn(1,(1+i)*tau)
                            res /= np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (1./np.pi**2)*quad(lambda o: np.exp( np.sqrt( t2+o**2 ) )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.)**2,0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def jp(self, tau: float, interp: bool= True) -> float: 
                """
                derivative j'(tau): cf. eq.(E.13)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series: 
                    zetapminus2 =-0.0304484570583933
                    return + tau*(7)*(zetapminus2)/np.pi**2
                else:
                    if interp:
                        return self.jp_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0.
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1.)**(i+1)*(i+1.)**2*kn(0,(1.+i)*tau)
                            res *= tau/np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (tau/np.pi**2)*quad(lambda o: np.exp( np.sqrt( t2+o**2 ) )\
                                    *( 1. - np.exp( np.sqrt( t2+o**2 ) ) )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.)**3/np.sqrt( t2+o**2 ),0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def J(self, tau: float, interp: bool= True) -> float:
                """
                cf. eq.(D.17) and eq.(E.14)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series: 
                    return 1/6. - tau**2/4/np.pi**2
                else:
                    if interp:
                        return self.J_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1)**i*kn(2,(1+i)*tau)
                            res *= tau**2/np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (1./np.pi**2)*quad(lambda o: o**2*np.exp( np.sqrt( t2+o**2 ) )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.)**2,0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def Y(self, tau: float, interp: bool= True) -> float:
                """
                cf. eq.(D.18) and eq.(E.15)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series: 
                    return 7*np.pi**2/30. - tau**2/4.
                else:
                    if interp:
                        return self.Y_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1)**i*kn(3,(1+i)*tau)/(i+1)
                            res *= 3*tau**3/np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (1./np.pi**2)*quad(lambda o: o**4*np.exp( np.sqrt( t2+o**2 ) )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.)**2,0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def k(self, tau: float, interp: bool= True) -> float:
                """
                cf. eq.(D.19) and eq.(E.16)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series: 
                    eulergamma  = 0.5772156649015329
                    return (np.log(np.pi/tau)-eulergamma)/(2*np.pi**2) 
                else:
                    if interp:
                        return self.k_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1)**i*kn(0,(1+i)*tau)
                            res /= np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (1./np.pi**2)*quad(lambda o: 1./np.sqrt( t2+o**2 )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.),0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def K(self, tau: float, interp: bool= True) -> float:
                """
                cf. eq.(D.20) and eq.(E.17)
                """
                if tau>NuDec_Const.tau_max_sf_series: 
                    return 0.
                elif tau<NuDec_Const.tau_min_sf_series/10:  # need an extra /10 to be safe 
                    return (1/4.)*( 1./3. + (tau/np.pi)**2*np.log(tau) )
                else:
                    if interp:
                        return self.K_interp(tau)
                    else:
                        if self.Bessel:
                            res = 0
                            for i in range(0,NuDec_Const.Bessel_Max):
                                res += (-1)**i*kn(1,(1+i)*tau)/(1+i)
                            res *= tau/np.pi**2
                            return res
                        else:
                            t2 = tau**2
                            return (1./np.pi**2)*quad(lambda o: o**2/np.sqrt( t2+o**2 )\
                                    /(np.exp( np.sqrt( t2+o**2 ) )+1.),0.,100.,epsabs=NuDec_Const.quad_atol,epsrel=NuDec_Const.quad_rtol)[0]


            def Z(self, tau: float, interp: bool= True) -> float:
                """
                cf. eqs.(D.21) and (D.22)
                """
                return 1/4*(self.Y(tau, interp) - 3*tau**2*self.K(tau, interp))



        class PolyLog:
            def __init__(self) -> None:
                pass

            def Li2(self, z: float) -> float:
                """
                polylogarithm function Li_2(z) using scipy.special
                """
                return spence(1.-z)

            def Li3(self, z: float) -> float:
                """
                polylogarithm function Li_3(z), for z = (-inf,1]
                """
                z3 = zeta(3)
                lz = np.log(abs(z))
                res = 0
                if (z<-1.):
                    res += self.Li3(1./z) - lz*( lz**2 + np.pi**2 )/6.
                elif (z==+1):
                    return z3
                elif (z==-1):
                    return (-3/4.)*z3
                elif (z<0):
                    return self.Li3(z**2)/4. - self.Li3(-z)
                elif (z<.25):
                    res = sum([ z**i/i**3 for i in range(1, 10)])
                elif (z>.25):
                    res  = z3 + lz*np.pi**2/6. + .5*lz**2*(1.5-np.log(-lz))
                    temp = 2.
                    for i in range(3,10):
                        temp *= i
                        res  += zeta(3-i)*lz**i/temp
                return res

            def Li4(self, z: float) -> float:
                """
                polylogarithm function Li_4(z), for z = (-inf,1]
                """
                z4 = zeta(4)
                lz = np.log(abs(z))
                res = 0
                if (z<-1.):
                    res -= self.Li4(1./z) + ( 7*np.pi**4 + 30*np.pi**2*lz**2 + 15*lz**4 )/360.
                elif (z==+1):
                    return z4
                elif (z==-1):
                    return (-7/8.)*z4
                elif (z<0):
                    return self.Li4(z**2)/8. - self.Li4(-z)
                elif (z<.25):
                    res = sum([ z**i/i**4 for i in range(1, 10)])
                elif (z>.25):
                    res  = z4 + (np.pi*lz)**2/12. + lz**3*(11/6.-np.log(-lz))/6.
                    res += lz*zeta(3)
                    temp = 6.
                    for i in range(4,10):
                        temp *= i
                        res  += zeta(4-i)*lz**i/temp
                return res



        # Pt_N: scaled coefficient functions: [like eq.(D.9), but without e^n part]

        def Pt_0(self, tau: float) -> float:
            return np.pi**2/45. + 2*self.sf.Z(tau,self.interp)/3.

        def Pt_2(self, tau: float) -> float:
            return -1./2.*( 1/3. + self.sf.K(tau,self.interp) )*self.sf.K(tau,self.interp)

        def Pt_3(self, tau: float) -> float:
            return self.sf.J(tau,self.interp)**(3/2)/3/np.sqrt(2.)/np.pi

        def Pt_4(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure: 
                return 0. # to avoid log(0)
            else:
                return -self.Pt_2(tau)*np.log(tau)/6/np.pi**2 + self.c4*np.pi**2/45.

        def Pt_5(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure:
                return 0. # to avoid log(0)
            else:
                return -self.Pt_3(tau)*np.log(tau)/4/np.pi**2 + self.c5*np.pi**2/45.


        # dPt_N: derivatives of Pt_N

        def dPt_0(self, tau: float) -> float:
            return -2*tau*self.sf.K(tau, self.interp)

        def dPt_2(self, tau: float) -> float:
            return tau*self.sf.k(tau,self.interp)*(1./6 + self.sf.K(tau,self.interp))

        def dPt_3(self, tau: float) -> float:
            return -tau*self.sf.j(tau,self.interp)*np.sqrt(self.sf.J(tau,self.interp))/2/np.sqrt(2.)/np.pi

        def dPt_4(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure: 
                return 0. # to avoid division by 0
            else:
                return ( self.sf.K(tau,self.interp)*(1.+3*self.sf.K(tau,self.interp))/tau \
                        - tau*self.sf.k(tau,self.interp)*np.log(tau)*(1.+6*self.sf.K(tau,self.interp)) )/36/np.pi**2

        def dPt_5(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure: 
                return 0. # to avoid log(0)
            else:
                return ( 3.*tau**2*self.sf.j(tau,self.interp)*np.log(tau) \
                    - 2.*self.sf.J(tau,self.interp) )*np.sqrt(self.sf.J(tau,self.interp))/24/np.sqrt(2)/np.pi**3/tau



        # G2_N: double derivatives of Pt_N [needed for d(Rho)/dT = 2*T^2*G_2(tau)]

        def G2_0(self, tau: float) -> float:
            return tau**2*self.sf.J(tau,self.interp) + self.sf.Y(tau,self.interp) + 2*np.pi**2/15.

        def G2_2(self, tau: float) -> float:
            return -( 3*self.sf.J(tau,self.interp) + 2*(self.sf.J(tau,self.interp)+self.sf.K(tau,self.interp))*(1.+3*self.sf.J(tau,self.interp)) \
                    + tau**2*self.sf.j(tau,self.interp)*(1.+6*self.sf.K(tau,self.interp)) )/12.

        def G2_3(self, tau: float) -> float:
            if tau>NuDec_Const.tau_max_sf_series: 
                j2_over_J   = 0
            else:
                j2_over_J = self.sf.j(tau,self.interp)**2/self.sf.J(tau,self.interp)

            return ( 8*self.sf.J(tau,self.interp) + 5*tau**2*self.sf.j(tau,self.interp) - tau**3*self.sf.jp(tau,self.interp) \
                + tau**4*j2_over_J/2. )*np.sqrt( self.sf.J(tau,self.interp)/2. )/np.pi/4.

        def G2_4(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure: 
                return 0. # to avoid log(0)
            else:
                return - self.G2_2(tau)*np.log(tau)/np.pi**2/6. + 2*self.c4*np.pi**2/15. \
                        - ( 3*self.sf.K(tau,self.interp)*(1-self.sf.K(tau,self.interp)) + 2*self.sf.J(tau,self.interp)*(1+6*self.sf.K(tau,self.interp)) )/72./np.pi**2

        def G2_5(self, tau: float) -> float:
            if tau<NuDec_Const.tau_min_QED_Pressure: 
                return 0. # to avoid log(0)
            else:
                return -self.G2_3(tau)*np.log(tau)/np.pi**2/4. + 2*self.c5*np.pi**2/15. \
                        + ( 7*self.sf.J(tau,self.interp) + 3*tau**2*self.sf.j(tau,self.interp) )*np.sqrt( self.sf.J(tau,self.interp)/2. )/np.pi**3/24.




        # bulk thermodynamic functions: [up to order e^5 in perturbation theory for QED]

        def P_QED_int(self, T: float, m: float=NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
            """
            perturbative QED pressure: cf. eq.(4.1)
            """

            if self.use_data_QED:
                funcs   = [self.p_int_2_non_ln, self.p_int_3, self.p_int_4, self.p_int_5]
                res     = self.p_int_0(T)
                res     += sum( funcs[i](T) for i in range(min(order-1, len(funcs))))
                if (order>5):
                    res     += self.p_int_2_ln(T) # separately include the log-term

                return res
            else:
                mT      = m / T
                res     = self.Pt_0(tau=mT)

                coeffs  = [e**2,e**3,e**4,e**5]
                funcs   = [self.Pt_2, self.Pt_3, self.Pt_4, self.Pt_5]

                res     += sum(coeffs[i] * funcs[i](tau=mT) for i in range(min(order-1, len(coeffs))))

                return res*T**4


        def I_QED_int(self, T: float, m: float=NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
            """
            perturbative `trace anomaly` rho-3P 
            """

            if self.use_data_QED:
                # I_QED_int = T*dp/dT - 4*p
                funcs   = [self.dp_int_2_non_ln, self.dp_int_3, self.dp_int_4, self.dp_int_5]
                res     = self.dp_int_0(T)
                res     += sum( funcs[i](T) for i in range(min(order-1, len(funcs))))
                if (order>5):
                    res     += self.dp_int_2_ln(T)  # separately include the log-term
                res     *= T
                res     -= 4*self.P_QED_int(T= T, m= NuDec_Const.me, order= order, e= NuDec_Const.e)

                return res
            else:
                mT      = m/T
                res     = self.dPt_0(tau=mT)

                coeffs  = [e**2,e**3,e**4,e**5]
                funcs   = [self.dPt_2, self.dPt_3, self.dPt_4, self.dPt_5]

                res     += sum(coeffs[i] * funcs[i](tau=mT) for i in range(min(order-1, len(coeffs))))

                return -m*res*T**3


        def Rho_QED_int(self, T: float, m: float=NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float: 
            """
            energy density: electrons & photons 
            """
            return self.I_QED_int(T=T,m=m,order=order,e=e)+3*self.P_QED_int(T=T,m=m,order=order,e=e)

        def s_QED_int(self, T: float, m: float, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float: 
            """
            entropy density: electrons & photons 
            """
            return ( self.I_QED_int(T=T, m=m, order=order, e=e)+4*self.P_QED_int(T=T, m=m, order=order, e=e) )/T


        def d2P_dT2_QED(self, T: float, m: float=NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float:
            """
            2nd derivative d^2P/dT^2 = (1/T)*dRho/dT 
            """

            if self.use_data_QED:
                funcs   = [self.d2p_int_2_non_ln, self.d2p_int_3, self.d2p_int_4, self.d2p_int_5]
                res     = self.d2p_int_0(T)
                res     += sum( funcs[i](T) for i in range(min(order-1, len(funcs))))
                if (order>5):
                    res     += self.d2p_int_2_ln(T) # separately include the log-term

                return res
            else:
                mT      = m/T
                res     = self.G2_0(tau=mT)

                coeffs  = [e**2,e**3,e**4,e**5]
                funcs   = [self.G2_2, self.G2_3, self.G2_4, self.G2_5]

                res     += sum(coeffs[i] * funcs[i](tau=mT) for i in range(min(order-1, len(coeffs))))

                return 2*res*T**2


        def P_FD(self, T: float, mu: float= 0, m: float= NuDec_Const.m_nu, g_internal: int= NuDec_Const.g_nu, Bessel: bool= True) -> float: 
            """
            pressure density of Fermi-Dirac species with g_internal degrees of freedom
            """
            lim_effectively_massless    = 1e-4
            if m/T<lim_effectively_massless:
                return -g_internal*T**4*self.polylog.Li4(-np.exp(mu/T))/np.pi**2
            else:
                if Bessel and mu/T<0.1:
                    n = np.arange(1, NuDec_Const.BesselMax + 1)  # vectorized n = 1..n_max
                    term = -g_internal*(((-1)**n) * np.exp(n * mu / T) * m**2 * T**2 * kv(2, m * n / T)) / (2 * n**2 * np.pi**2)
                    return np.sum(term)
                else:
                    if Bessel and mu/T>=0.1:
                        warnings.warn("Bessel expansion requested in P_FD, but mu/T>0.1 for it is not convergent. " \
                                        "Fall back to numerical integration..", category=UserWarning, stacklevel=2)
                    return g_internal*quad(lambda k: 1/(6*np.pi**2)*k**4/np.sqrt(k**2+m**2)*1/(np.exp( np.sqrt(k**2+m**2)/T - mu/T) + 1), 0, 25*T+mu)[0]



        def n_FD(self, T: float, mu: float= 0, m: float= NuDec_Const.m_nu, g_internal: int= NuDec_Const.g_nu, Bessel: bool= True) -> float: 
            """
            number density of Fermi-Dirac species with g_internal degrees of freedom
            """
            lim_effectively_massless    = 1e-4
            if m/T<lim_effectively_massless:
                return -g_internal*T**3*self.polylog.Li3(-np.exp(mu/T))/np.pi**2 
            else:
                if Bessel and mu/T<0.1:
                    n       = np.arange(1, NuDec_Const.BesselMax + 1)  
                    term    = -g_internal*(((-1)**n) * np.exp(n * mu / T) * m**2 * T * kv(2, m * n / T)) / (2 * n * np.pi**2)
                    return np.sum(term)
                else:
                    if Bessel and mu/T>=0.1:
                        warnings.warn("Bessel expansion requested in n_FD, but mu/T>0.1 for it is not convergent. " \
                                        "Fall back to numerical integration..", category=UserWarning, stacklevel=2)
                    return g_internal*quad(lambda k: 1/(2*np.pi**2)*k**2*1/(np.exp( np.sqrt(k**2+m**2)/T - mu/T) + 1), 0, 25*T+mu)[0]




        def Rho_FD(self, T: float, mu: float= 0, m: float= NuDec_Const.m_nu, g_internal: int= NuDec_Const.g_nu, Bessel: bool= True) -> float: 
            """
            pressure density of Fermi-Dirac species with g_internal degrees of freedom
            """
            lim_effectively_massless    = 1e-4
            if m/T<lim_effectively_massless:
                return 3.*self.P_FD(T,mu,m,g_internal,Bessel)
            else:
                if Bessel and mu/T<0.1:
                    n = np.arange(1, NuDec_Const.BesselMax + 1)
                    x = m * n / T
                    term = -g_internal*(((-1)**n) * np.exp(n * mu / T) * m**2 * T *
                            (m * n * kv(1, x) + 3 * T * kv(2, x))) / (2 * n**2 * np.pi**2)
                    return np.sum(term)
                else:
                    if Bessel and mu/T>=0.1:
                        warnings.warn("Bessel expansion requested in Rho_FD, but mu/T>0.1 for it is not convergent. " \
                                        "Fall back to numerical integration..", category=UserWarning, stacklevel=2)
                    return g_internal*quad(lambda k: 1/(2*np.pi**2)*k**2*np.sqrt(k**2+m**2)*1/(np.exp( np.sqrt(k**2+m**2)/T - mu/T) + 1), 0, 25*T+mu)[0]



        def s_FD(self, T: float, mu: float= 0, m: float= NuDec_Const.m_nu, g_internal: int= NuDec_Const.g_nu, Bessel: bool= True) -> float: 
            """
            entropy density of neutrinos
            """
            return ( self.P_FD(T=T, mu=mu, m=m, g_internal= g_internal, Bessel= Bessel) + self.Rho_FD(T=T, mu=mu, m=m, g_internal= g_internal, Bessel= Bessel) - mu*self.n_FD(T=T, mu=mu, m=m, g_internal= g_internal, Bessel= Bessel) )/T


        # thermodynamics for neutrinos, cf. table V
        def dn_dT_nu(self, T: float, mu: float= 0) -> float:
            """
            1st derivative dn/dT
            """
            return -2.*T*( 3*T*self.polylog.Li3(-np.exp(mu/T)) - mu*self.polylog.Li2(-np.exp(mu/T)) )/np.pi**2

        def dn_dmu_nu(self, T: float, mu: float= 0) -> float:
            """
            1st derivative dn/dmu
            """
            return -2.*(T/np.pi)**2*self.polylog.Li2(-np.exp(mu/T))

        def dP_dT_nu(self, T: float, mu: float= 0) -> float:
            """
            1st derivative dP/dT
            """
            return -2.*T**2*( 4*T*self.polylog.Li4(-np.exp(mu/T)) - mu*self.polylog.Li3(-np.exp(mu/T)) )/np.pi**2

        def dP_dmu_nu(self, T: float, mu: float= 0) -> float:
            """
            1st derivative dP/dmu
            """
            return -2.*T**3*self.polylog.Li3(-np.exp(mu/T))/np.pi**2

        def Rho_tot(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float: 
            """
            total energy density
            """
            return self.Rho_QED_int(T=T_gam, m=me, order=order, e=e) + self.Rho_FD(T= T_nue, mu= mu_nue, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu) + 2*self.Rho_FD(T=T_numu, mu= mu_numu, m=NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu)

        def P_tot(self, T_gam: float, T_nue: float, T_numu: float, mu_nue: float, mu_numu: float, me: float= NuDec_Const.me, order: int= NuDec_Const.order, e: float= NuDec_Const.e) -> float: 
            """
            total pressure 
            """
            return self.P_QED_int(T=T_gam, m= me, order=order, e=e) + self.P_FD(T=T_nue, mu=mu_nue, m=NuDec_Const.m_nu, g_internal= NuDec_Const.g_nu) + 2*self.P_FD(T=T_numu, mu=mu_numu, m=NuDec_Const.m_nu, g_internal=NuDec_Const.g_nu)


