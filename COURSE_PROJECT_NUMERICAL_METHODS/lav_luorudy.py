# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 17:12:51 2023

@author: Kaviyarasan PR
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
Cm = 1.0  # Membrane capacitance (uF/cm^2)
ENa = 54.4  # Sodium equilibrium potential (mV)
EK = -77  # Potassium equilibrium potential (mV)
EL = -54.4  # Leak equilibrium potential (mV)
GNa = 23  # Sodium conductance (mS/cm^2)
GK = 0.282  # Potassium conductance (mS/cm^2)
GL = 0.3  # Leak conductance (mS/cm^2)

# Additional constants for potassium currents
GK1 = 0.6047 * (np.sqrt(5.4) / 5.4)  # Time-independent potassium conductance (mS/cm^2)
GKP = 0.282 * (np.sqrt(5.4) / 5.4)  # Plateau potassium conductance (mS/cm^2)
Gb = 0.0392  # Background potassium conductance (mS/cm^2)
EK1 = -87.2  # Time-independent potassium equilibrium potential (mV)
EKp = -87.2  # Plateau potassium equilibrium potential (mV)
Eb = -59.87  # Background potassium equilibrium potential (mV)

# Constants for slow inward current
Gsi = 0.09  # Slow inward current conductance (mS/cm^2)
Esi = 118.7  # Slow inward current equilibrium potential (mV)

# Simulation parameters
dt = 0.01  # Time step (ms)
tmax = 500  # Maximum time (ms) - Increased to 500
t = np.arange(0, tmax, dt)  # Time vector

# Initial conditions
V = -84  # Membrane potential (mV)
n = 0.317  # Potassium gating variable n
m = 0.05  # Sodium activation gating variable m
h = 0.6  # Sodium inactivation gating variable h
X = 0.1  # Activation gate for time-dependent potassium current
Xi = 0.1  # Inactivation gate for time-dependent potassium current
K1_inf = 0.1  # Inactivation gate for time-independent potassium current
Kp = 0.1  # Plateau potassium current gating variable
d = 0.0  # Activation gate for slow inward current
f = 1.0  # Inactivation gate for slow inward current
Cai = 0.0001  # Initial calcium concentration (mM)
j= 0.1

# Lists to store results
V_trace = [V]
n_trace = [n]
m_trace = [m]
h_trace = [h]
X_trace = [X]
Xi_trace = [Xi]
K1_inf_trace = [K1_inf]
Kp_trace = [Kp]
j_trace = [j]
d_trace = [d]
f_trace = [f]
Cai_trace = [Cai]

# Simulation
for i in range(1, len(t)):
    # Compute gating variable changes
    alpha_m = (0.32 * (V + 47.13)) / (1 - np.exp(-(V + 47.13) / 10))
    beta_m = 0.08 * np.exp(-(V) / 11)
    alpha_h = 0.135 * np.exp((V + 80) / -6.8) if V < -40 else 0
    beta_h = (3.56 * np.exp(0.079 * V) + 3.1e5 * np.exp(0.35 * V)) if V < -40 else (1 / (0.13 * (1 + np.exp((V + 10.66) / -11.1))))
    alpha_n = (0.01 * (V + 55)) / (1 - np.exp(-(V + 55) / 10))
    beta_n = 0.125 * np.exp(-(V + 65) / 80)
    alpha_j = ((-1.2714e5 * np.exp(0.2444 * V)) - (3.474e-5 * np.exp(-0.04391 * V)) * (V + 37.78))/(1+ np.exp(0.311*(V + 79.23))) if V < -40 else 0
    beta_j = (0.1212 * np.exp(-0.01052 * V)) / (1 + np.exp(-0.1378 * (V + 40.14))) if V < -40 else (0.3 * np.exp(-2.535e-7 * V)) / (1 + np.exp(-0.1 * (V + 32)))

    # Calculate time-dependent potassium current (IK)
    alpha_X = 0.0005 * np.exp(0.083 * (V + 50)) / (1 + np.exp(0.057 * (V + 50)))
    beta_X = 0.0013 * np.exp(-0.06 * (V + 20)) / (1 + np.exp(-0.04 * (V + 20)))
    
    # Compute changes for the slow inward current gates (d and f)
    alpha_d = 0.095 * np.exp(-0.01 * (V - 5)) / (1 + np.exp(-0.072 * (V - 5)))
    beta_d = 0.07 * np.exp(-0.017 * (V + 44)) / (1 + np.exp(0.05 * (V + 44)))
    alpha_f = 0.012 * np.exp(-0.008 * (V + 28)) / (1 + np.exp(0.15 * (V + 28)))
    beta_f = 0.0065 * np.exp(-0.02 * (V + 30)) / (1 + np.exp(-0.2 * (V + 30)))

    dd = alpha_d * (1 - Cai) - beta_d * Cai
    df = alpha_f * (1 - Cai) - beta_f * Cai

   

    # Add the equations for Xi, K1_inf, and Kp here
    if V > -100:
        Xi = (2.837 * (np.exp(0.04 * (V + 77)) - 1) * (V + 77) * np.exp(0.04 * (V + 35))) / (1 + np.exp(0.04 * (V + 35)))
    else:
        Xi = 1
    alpha_K1 = 1.02 / (1 + np.exp(0.2385 * (V - EK1 - 59.215)))
    beta_K1 = (0.49124 * np.exp(0.08032 * (V - EK1 + 5.476)) + np.exp(0.06175 * (V - EK1 - 594.31))) / (1 + np.exp(-0.5143 * (V - EK1 + 4.753)))
    K1_inf = alpha_K1 / (alpha_K1 + beta_K1)
    Kp = 1 / (1 + np.exp((7.488 - V) / 5.98))

    # Compute total ionic current Iion
    

    # Update Ca_i based on the calcium uptake equation
    dCai = -1e-4 * Isi + 0.07 * (1e-4 - Cai)
    Cai += dCai * dt
    
    
    
    INa = GNa * m**3 * h * (V - ENa)
    IK = GK * n**4 * (V - EK)
    IL = GL * (V - EL)
    IK1 = GK1 * K1_inf * (V - EK1)
    IKp = GKP * Kp * (V - EKp)
    Ib = Gb * (V - Eb)
    Isi = Gsi * d * f * (V - Esi)
    Iion = INa + IK + IL + IK1 + IKp + Ib + Isi

    # Update membrane potential and gating variables
    dV = (1 / Cm) * (-(Iion+ 1.00005555))
    V += dV * dt
    dm = (0.32 * (V + 47.13)) / (1 - np.exp(-(V + 47.13) / 10)) * (1 - m) - 0.08 * np.exp(-(V) / 11) * m
    dj = ((-1.2714e5 * np.exp(0.2444 * V)) - (3.474e-5 * np.exp(-0.04391 * V)) * (V + 37.78))/(1+ np.exp(0.311*(V + 79.23))) * (1 - j) - (0.1212 * np.exp(-0.01052 * V)) / (1 + np.exp(-0.1378 * (V + 40.14))) * j
    dh = 0.135 * np.exp((V + 80) / -6.8) * (1 - h) - (3.56 * np.exp(0.079 * V) + 3.1e5 * np.exp(0.35 * V)) * h
    dn = (0.01 * (V + 55)) / (1 - np.exp(-(V + 55) / 10)) * (1 - n) - 0.125 * np.exp(-(V + 65) / 80) * n
    dX = 0.0005 * np.exp(0.083 * (V + 50)) / (1 + np.exp(0.057 * (V + 50))) * (1 - X) - 0.0013 * np.exp(-0.06 * (V + 20)) / (1 + np.exp(-0.04 * (V + 20))) * X
    dd = 0.095 * np.exp(-0.01 * (V - 5)) / (1 + np.exp(-0.072 * (V - 5))) * (1 - d) - 0.07 * np.exp(-0.017 * (V + 44)) / (1 + np.exp(0.05 * (V + 44)))* d
    df = 0.012 * np.exp(-0.008 * (V + 28)) / (1 + np.exp(0.15 * (V + 28))) * (1 - f) - 0.0065 * np.exp(-0.02 * (V + 30)) / (1 + np.exp(-0.2 * (V + 30))) * f

    n += dn * dt
    m += dm * dt
    h += dh * dt
    j += dj * dt
    X += dX * dt
    d += dd * dt
    f += df * dt

    # Append results to trace lists
    V_trace.append(V)
    n_trace.append(n)
    m_trace.append(m)
    h_trace.append(h)
    X_trace.append(X)
    Xi_trace.append(Xi)
    K1_inf_trace.append(K1_inf)
    Kp_trace.append(Kp)
    d_trace.append(d)
    f_trace.append(f)
    Cai_trace.append(Cai)

# Plot the transmembrane potential with y-axis ranging from -100 to 50
plt.figure(figsize=(12, 8))
plt.plot(t, V_trace, 'b', linewidth=2)
plt.title('Transmembrane Potential (Action Potential)')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.ylim(-100, 50)  # Set the y-axis limits
plt.grid()
plt.tight_layout()
plt.show()
