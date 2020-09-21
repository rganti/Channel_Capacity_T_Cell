from pysb import *

from simulation_parameters import BindingParameters

rate_constants = BindingParameters()

Model()

Monomer('R', ['b', 's'], {'s': ['u', 'p']})
Monomer('Ls', ['b'])
Monomer('Lf', ['b'])

Parameter('k_on_L', rate_constants.k_L_on)
Parameter('k_off_Lf', rate_constants.k_foreign_off)
Parameter('k_off_Ls', rate_constants.k_self_off)

Parameter('R_0', 10000)
Parameter('Ls_0', 800)
Parameter('Lf_0', 30)

Initial(R(b=None, s='u'), R_0)
Initial(Lf(b=None), Lf_0)
Initial(Ls(b=None), Ls_0)

Rule('R_Lf_bind', R(b=None) + Lf(b=None) | R(b=1) % Lf(b=1), k_on_L, k_off_Lf)
Observable("O_Lf", Lf(b=None))
Observable("O_RLf", R(b=1) % Lf(b=1))

Rule('R_Ls_bind', R(b=None) + Ls(b=None) | R(b=1) % Ls(b=1), k_on_L, k_off_Ls)
Observable("O_Ls", Ls(b=None))
Observable("O_RLs", R(b=1) % Ls(b=1))
