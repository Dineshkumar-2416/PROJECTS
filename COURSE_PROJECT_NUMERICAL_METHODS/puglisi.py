# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 18:21:44 2023

@author: Kaviyarasan PR
"""
INa= GNa*(m**3)*h*j*(V-ENa)
import numpy as np
from scipy.integrate import odeint

def cardiac_cell_model(y, t, Istim, A_cap, V_myo, F, V_JSR, V_NSR, RAV):
    Vm, m, h, j, d, f, b, g, Xr, Xs, Na_i, Ca_i, K_i, Ca_JSR, Ca_NSR, Ca_foot = y

    # Parameters (You should provide appropriate values for these)
    # Define all the parameters, such as alpha, beta, and other constants here.

    # Differential equations
    dVm_dt = (Istim - (INa + ICaL + ICaT + IKr + IKs + INaCa + IK1 + IKp) / C + Istim - (IpCa + INab + ICab + INaK + Ito + ICl_Ca) / C) / 2

    dm_dt = alpha_m * (1 - m) - beta_m * m
    dh_dt = alpha_h * (1 - h) - beta_h * h
    dj_dt = alpha_j * (1 - j) - beta_j * j
    dd_dt = alpha_d * (1 - d) - beta_d * d
    df_dt = alpha_f * (1 - f) - beta_f * f
    db_dt = (b_inf - b) / tau_b
    dg_dt = (g_inf - g) / tau_g
    dXr_dt = (Xr_inf - Xr) / tau_Xr
    dXs_dt = (Xs_inf - Xs) / tau_Xs

    dNa_i_dt = - (INa + ICaNa + INab + 3 * INaCa + 3 * INaK) * A_cap / (V_myo * F)
    dCa_i_dt = ((ICaCa + IpCa + ICab + ICaT) - INaCa) * A_cap / (2 * V_myo * F) + Irel * V_JSR / V_myo + (Ileak - Iup) * V_NSR / V_myo
    dK_i_dt = - (ICaK + IKr + IKs + IK1 + IKp + Ito - 2 * INaK) * A_cap / (V_myo * F)

    dCa_JSR_dt = - (Irel - Itr * V_NSR / V_JSR)
    dCa_NSR_dt = - ((Ileak + Itr) - Iup)
    dCa_foot_dt = (ICaCa) * (A_cap / (2 * V_myo * F)) * RAV

    dydt = [dVm_dt, dm_dt, dh_dt, dj_dt, dd_dt, df_dt, db_dt, dg_dt, dXr_dt, dXs_dt, dNa_i_dt, dCa_i_dt, dK_i_dt, dCa_JSR_dt, dCa_NSR_dt, dCa_foot_dt]
    return dydt

# Example usage
# Define initial conditions and other parameters
initial_conditions = [Vm0, m0, h0, j0, d0, f0, b0, g0, Xr0, Xs0, Na_i0, Ca_i0, K_i0, Ca_JSR0, Ca_NSR0, Ca_foot0]
time = np.linspace(0, 100, 1000)  # Time vector
Istim = 1.0  # Example stimulus current

# Solve the differential equations
solution = odeint(cardiac_cell_model, initial_conditions, time, args=(Istim, A_cap, V_myo, F, V_JSR, V_NSR, RAV))

# Access the results
Vm, m, h, j, d, f, b, g, Xr, Xs, Na_i, Ca_i, K_i, Ca_JSR, Ca_NSR, Ca_foot = solution.T
