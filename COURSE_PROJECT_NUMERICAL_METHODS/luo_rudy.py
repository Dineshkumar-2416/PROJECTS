
###trie///////LUO RUDY ELECTRO PHYSIOLOGICAL MODEL FOR CARDIO /////////

from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.misc import derivative
from operator import itemgetter
from scipy.optimize import newton
from matplotlib.ticker import MultipleLocator, FuncFormatter
from sympy import Eq, dsolve, Function
from scipy.integrate import ode



# In[2]:


def func(w, t):
    V, m, h, n, X, d, f = w
    dm = lambda m, t, V:  (0.32 * (V + 47.13)) / (1 - np.exp(-(V + 47.13) / 10)) * (1 - m) - 0.08 * np.exp(-(V) / 11) * m
    dh = lambda h, t, V: 0.135 * np.exp((V + 80) / -6.8) * (1 - h) - (3.56 * np.exp(0.079 * V) + 3.1e5 * np.exp(0.35 * V)) * h
    dn = lambda n, t, V: (0.01 * (V + 55)) / (1 - np.exp(-(V + 55) / 10)) * (1 - n) - 0.125 * np.exp(-(V + 65) / 80) * n
    dX=  lambda X ,t, V: 0.0005 * np.exp(0.083 * (V + 50)) / (1 + np.exp(0.057 * (V + 50))) * (1 - X) - 0.0013 * np.exp(-0.06 * (V + 20)) / (1 + np.exp(-0.04 * (V + 20))) * X
    dd= lambda d, t, V: 0.095 * np.exp(-0.01 * (V - 5)) / (1 + np.exp(-0.072 * (V - 5))) * (1 - d) - 0.07 * np.exp(-0.017 * (V + 44)) / (1 + np.exp(0.05 * (V + 44)))* d
    df= lambda f ,t, V: 0.012 * np.exp(-0.008 * (V + 28)) / (1 + np.exp(0.15 * (V + 28))) * (1 - f) - 0.0065 * np.exp(-0.02 * (V + 30)) / (1 + np.exp(-0.2 * (V + 30))) * f
    dV = lambda V, t, m, h, n, X, d, f : -((GNa * m**3 * h * (V - ENa)+ GK * n**4 * (V - EK) + GL * (V - EL) + GK1 * K1_inf * (V - EK1) + GKP * Kp * (V - EKp) + Gb * (V - Eb) + Gsi * d * f * (V - Esi))+Ist)/12
    return np.array([dV(V,t,m,h,n,X,d,f),dm(m,t,V),dh(h,t,V),dn(n,t,V),dX(X,t,V),dd(d,t,V),df(f,t,V)])


# ## Kinetics kurve

def runge_kutta_4(func, w0, t):
    # Initialize arrays to store the solutions
    num_points = len(t)
    w = np.zeros((num_points, len(w0)))
    w[0] = w0

    # Perform RK4 integration
    for i in range(num_points - 1):
        h = t[i + 1] - t[i]
        k1 = h * func(w[i], t[i])
        k2 = h * func(w[i] + 0.5 * k1, t[i] + 0.5 * h)
        k3 = h * func(w[i] + 0.5 * k2, t[i] + 0.5 * h)
        k4 = h * func(w[i] + k3, t[i + 1])
        w[i + 1] = w[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return w

def implicit_rk4(func, y0, t):
    n = len(t)
    m = len(y0)
    y = np.zeros((n, m))
    y[0] = y0

    for i in range(n - 1):
        h = t[i + 1] - t[i]
        k1 = func(y[i], t[i])
        
        # Define a function for Newton's method to find the new state
        def equation(y_new):
            k2 = func(y_new, t[i] + h / 2)
            k3 = func(y_new, t[i] + h / 2)
            k4 = func(y_new, t[i] + h)
            y_next = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
            return y[i + 1] - y_new + y_next
        
        # Use Newton's method to find the new state
        y_new = newton(equation, y[i])
        y[i + 1] = y_new

    return y

def euler(func, w0, t):
    # Initialize arrays to store the solutions
    num_points = len(t)
    w = np.zeros((num_points, len(w0)))
    w[0] = w0

    # Perform RK4 integration
    for i in range(num_points - 1):
        h = t[i + 1] - t[i]
        k1 = h * func(w[i], t[i])
        w[i + 1] = w[i] + k1

    return w


# In[3]:


Cm = 1.0  # Membrane capacitance (uF/cm^2)
ENa = 54.4  # Sodium equilibrium potential (mV)
EK = -77  # Potassium equilibrium potential (mV)
EL = -54.4  # Leak equilibrium potential (mV)
GNa = 23  # Sodium conductance (mS/cm^2)
GK = 0.282  # Potassium conductance (mS/cm^2)
GL = 0.3  # Leak conductance (mS/cm^2)

GK1 = 0.6047 * (np.sqrt(5.4) / 5.4)  # Time-independent potassium conductance (mS/cm^2)
GKP = 0.282 * (np.sqrt(5.4) / 5.4)  # Plateau potassium conductance (mS/cm^2)
Gb = 0.0392  # Background potassium conductance (mS/cm^2)
EK1 = -87.2  # Time-independent potassium equilibrium potential (mV)
EKp = -87.2  # Plateau potassium equilibrium potential (mV)
Eb = -59.87  # Background potassium equilibrium potential (mV)


Gsi = 0.09  # Slow inward current conductance (mS/cm^2)
Esi = 118.7  # Slow inward current equilibrium potential (mV)


K1_inf = 0.1  # Inactivation gate for time-independent potassium current
Kp = 0.1  # Plateau potassium current gating variable
Ist= 5 



t = np.linspace(0, 1000, 100000)
sol = euler(func, (-90,0.07,0.06,0.03,0.1,0.0,1.0), t)
sol1 = implicit_rk4(func, (-80,0.07,0.06,0.03,0.1,0.0,1.0), t)
sol2 = odeint(func, (-80,0.07,0.06,0.03,0.1,0.0,1.0), t)
sol3 = runge_kutta_4(func, (-80,0.07,0.06,0.03,0.1,0.0,1.0), t)


# ... (Previous code remains the same)

# Plot the transmembrane potential with higher precision in the y-axis
plt.figure(figsize=(12, 8))
plt.plot(t, sol[:,0], linewidth=2 ,label= 'Euler')
plt.plot(t,sol1[:,0],label ='Implicit rk4')
plt.plot(t,sol2[:,0],label ='Inbuilt odeint')
plt.plot(t,sol3[:,0],label ='Runge kutta 4')

plt.title('Transmembrane Potential (Action Potential)')

plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')


plt.grid()

# Customize the y-axis to show more decimal places
plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.7f'))
plt.gca().yaxis.set_major_locator(plt.MultipleLocator(1))  # You can adjust the increment as needed

plt.tight_layout()
plt.show()