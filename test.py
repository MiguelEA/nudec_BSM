import time
import numpy as np
import matplotlib.pyplot as plt
from importlib import reload

import source.nudec_source as ns

nudec = ns.NuDec()

start_time = time.time()

[t, y]     = nudec.evolve()
# output array: y = [T_gam, T_nue, T_numu, mu_nue, mu_numu, z]

end_time = time.time()
print("")
print("solver time:", f"{(end_time - start_time):.3g}", "[s]")

# access final entries in output:
_Tg = y[0][-1]
_Te = y[1][-1]
_Tm = y[2][-1]
_me = y[3][-1]
_mm = y[4][-1]
_z  = y[5][-1]

# compute observables:
neff         =         nudec.Neff(_Tg, _Te, _Tm, _me, _mm)
gstar_rho    =    nudec.gstar_rho(_Tg, _Te, _Tm, _me, _mm)
gstar_s      =      nudec.gstar_s(_Tg, _Te, _Tm, _me, _mm)
m_over_Omega = nudec.m_over_Omega(_Tg, _Te, _Tm, _me, _mm)
heff         = nudec.heff(_z)

print("Neff         = ", f"{neff:.10g}")
print("gstar_rho    = ", f"{gstar_rho:.10g}")
print("gstar_s      = ", f"{gstar_s:.10g}")
print("m_over_Omega = ", f"{m_over_Omega:.10g}")
print("heff         = ", f"{heff:.10g}")


# write thermodynamic history to file:A
fname = "scan_hist.dat"
with open(fname, "w") as fout:
    fout.write("# Columns: t (s), T_gam (MeV), T_nue (MeV), T_numu (MeV), mu_nue (MeV), mu_numu (MeV), z\n")

    for i in range(len(t)):

        _t  = t[i]
        _Tg = y[0][i]
        _Te = y[1][i]
        _Tm = y[2][i]
        _me = y[3][i]
        _mm = y[4][i]
        _z  = y[5][i]

        fout.write(f"{_t: .8e}   {_Tg: .8e}   {_Te: .8e}   {_Tm: .8e}   {_me: .8e}  {_mm: .8e}   {_z: .8e}\n")

    print("\n --> output written to file: ["+fname+"]\n")

