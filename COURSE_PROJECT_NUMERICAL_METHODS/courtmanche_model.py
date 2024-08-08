# tried///////COURTMANCHE ELECTRO PHYSIOLOGICAL MODEL/////



import numpy as np

# Define the constants and parameters as global variables
R = 8.3143  # J·K^-1·mol^-1
T = 310  # K
F = 96.4867  # C/mmol
Cm = 100  # pF
Vcell = 20100  # μm^3
Vi = 13668  # μm^3
Vup = 1109.52  # μm^3
Vrel = 96.48  # μm^3
K1o = 5.4  # mM
Na1o = 140  # mM
Ca21o = 1.8  # mM
gNa = 7.8  # nS/pF
gK1 = 0.09  # nS/pF
gto = 0.1652  # nS/pF
gKr = 0.0294  # nS/pF
gKs = 0.129  # nS/pF
gCaL = 0.1238  # nS/pF
gbCa = 0.00113  # nS/pF
gbNa = 0.000674  # nS/pF
INaK_max = 0.60  # pA/pF
INaCa_max = 1600  # pA/pF
IpCa_max = 0.275  # pA/pF
Iup_max = 0.005  # mM/ms
KQ10 = 3
g = 0.35
Km_Nai = 10  # mM
Km_Ko = 1.5  # mM
Km_Na = 87.5  # mM
Km_Ca = 1.38
ksat = 0.1
krel = 30e-3  # 30 ms^-1
Kup = 0.00092  # mM
Ca21_up_max = 15  # mM
Cmdn_max = 0.05  # mM
Trpn_max = 0.07  # mM
Csqn_max = 10  # mM
Km_Cmdn = 0.00238  # mM
Km_Trpn = 0.0005  # mM
Km_Csqn = 0.8  # mM

def model(w, t):
    global R, T, F, Cm, Vcell, Vi, Vup, Vrel, K1o, Na1o, Ca21o, gNa, gK1, gto, gKr, gKs, gCaL, gbCa, gbNa, INaK_max, INaCa_max, IpCa_max, Iup_max, KQ10, g, Km_Nai, Km_Ko, Km_Na, Km_Ca, ksat, krel, Kup, Ca21_up_max, Cmdn_max, Trpn_max, Csqn_max, Km_Cmdn, Km_Trpn, Km_Csqn

    V, m, h, j, oa, oi, ua, ui, xr, xs, d, f, f_Ca, u, v, w, Ca_i, Ca_rel, Ca_up, Na_i, K_i = w

    # Constants
    RTF = R * T / F
    E_Na = RTF * np.log(Na1o / Na_i)
    E_K = RTF * np.log(K1o / K_i)
    E_Ca = RTF * np.log(Ca21o / Ca_i)
    sigma = (1/7) * (np.exp(Na1o / 67.3) - 1)
    Fn = (10**(-12))*Vrel*i_rel - (5*(10**(-13)))*(i_Ca_L/2 - i_NaCa/5 )/F

    # Component: fast sodium current (i_Na)
    i_Na = Cm * g_Na * m**3 * h * j * (V - E_Na)

    # Component: fast sodium current m gate
    alpha_m = np.where(np.abs(V + 47.13) < 1e-6, 3.2, 0.32 * (V + 47.13) / (1 - np.exp(-0.1 * (V + 47.13)))

    beta_m = 0.08 * np.exp(-V / 11)
    m_inf = alpha_m / (alpha_m + beta_m)
    tau_m = 1 / (alpha_m + beta_m)
    dmdt = (m_inf - m) / tau_m

    # Component: fast sodium current h gate
    alpha_h = np.where(V < -400, 0.135 * np.exp(V + 80) - 6.8, 3.56 * np.exp(0.079 * V + 3.1) / (5 * np.exp(0.35 * V)))
    beta_h = np.where(V < -400, 0.3 * np.exp(-2.535e-7 * V) / (1 + np.exp(-0.1 * (V + 32))),
                     0.13 / (1 + np.exp(V + 10.66 - 11.1)))
    h_inf = alpha_h / (alpha_h + beta_h)
    tau_h = 1 / (alpha_h + beta_h)
    dhdt = (h_inf - h) / tau_h

    # Component: fast sodium current j gate
    alpha_j = np.where(V < -400, ((-1.2714e5 * np.exp(0.2444 * V) - 3.474e-5 * np.exp(-0.04391 * V)) * (V + 37.78) / (1 + np.exp(0.311 * (V + 79.23))),
                     0.1212 * np.exp(-0.01052 * V) / (1 + np.exp(-0.1378 * (V + 40.14)))
    beta_j = np.where(V < -400, 0.3 * np.exp(-2.535e-7 * V) / (1 + np.exp(-0.1 * (V + 32))),
                    0.0  # Modify the else part accordingly
    j_inf = alpha_j / (alpha_j + beta_j)
    tau_j = 1 / (alpha_j + beta_j)
    djdt = (j_inf - j) / tau_j

    # Component: transient outward K current (i_to)
    i_to = Cm * g_to * oa * (oi + 0.5 * u) * (V - E_K)

    # Component: transient outward K current oa gate
    alpha_oa = 0.65 * (np.exp(V - (-10)) - 8.5 + np.exp(V - (-10)) - 40 - 59)**-1
    beta_oa = 0.65 * (2.5 + np.exp(V - (-10)) + 7217)**-1
    tau_oa = (alpha_oa + beta_oa)**-1 * K_Q10
    oa_infinity = (1 + np.exp(V - (-10)) + 10.47 - 17.54)**-1
    doa_dt = (oa_infinity - oa) / tau_oa

    # Component: transient outward K current oi gate
    alpha_oi = (18.53 + 1 * np.exp(V - (-10)) + 103.7) / (10.95)**-1
    beta_oi = (35.56 + 1 * np.exp(V - (-10)) - 8.74 - 7.44)**-1
    tau_oi = (alpha_oi + beta_oi)**-1 * K_Q10
    oi_infinity = (1 + np.exp(V - (-10)) + 33.15 / 5.3)**-1
    doi_dt = (oi_infinity - oi) / tau_oi

    # Component: ultrarapid delayed rectifier K current (i_Kur)
    g_Kur = 0.005 + 0.051 / (np.exp(V - 15) - 13)
    i_Kur = Cm * g_Kur * ua * (ui + 0.5 * u) * (V - E_K)

    # Component: ultrarapid delayed rectifier K current ua gate
    alpha_ua = 0.65 * (np.exp(V - (-10)) - 8.5 + np.exp(V - (-10)) - 40 - 59)**-1
    beta_ua = 0.65 * (2.5 + np.exp(V - (-10)) + 7217)**-1
    tau_ua = (alpha_ua + beta_ua)**-1 * K_Q10
    ua_infinity = (1 + np.exp(V - (-10)) + 20.3 - 9.6)**-1
    dua_dt = (ua_infinity - ua) / tau_ua

    # Component: ultrarapid delayed rectifier K current ui gate
    alpha_ui = (21 + 1 * np.exp(V - (-10)) - 195 - 28)**-1
    beta_ui = 1 * np.exp(V - (-10)) - 168 - 16
    tau_ui = (alpha_ui + beta_ui)**-1 * K_Q10
    ui_infinity = (1 + np.exp(V - (-10)) - 109.45 / 27.48)**-1
    dui_dt = (ui_infinity - ui) / tau_ui

    # Component: rapid delayed rectifier K current (i_Kr)
    i_Kr = Cm * g_Kr * xr * (V - E_K) / (1 + np.exp(V + 15) / 22.4)

    # Component: rapid delayed rectifier K current xr gate
    alpha_xr = np.where(abs(V + 14.1) < 1e-10, 0.0015,
                       0.0003 * (V + 14.1) / (1 - np.exp(V + 14.1 - 5))
    beta_xr = np.where(abs(V - 3.3328) < 1e-10, 3.7836118e-4,
                      0.000073898 * (V - 3.3328) * np.exp(V - 3.3328) / (5.1237 - 1)
    tau_xr = (alpha_xr + beta_xr)**-1
    xr_infinity = (1 + np.exp(V + 14.1 - 6.5))**-1
    dxr_dt = (xr_infinity - xr) / tau_xr

    # Component: slow delayed rectifier K current (i_Ks)
    i_Ks = Cm * g_Ks * xs**2 * (V - E_K)

    # Component: slow delayed rectifier K current xs gate
    alpha_xs = np.where(abs(V - 19.9) < 1e-10, 0.00068,
                       0.00004 * (V - 19.9) / (1 - np.exp(V - 19.9 - 17))
    beta_xs = np.where(abs(V - 19.9) < 1e-10, 0.000315,
                      0.000035 * (V - 19.9) * np.exp(V - 19.9) / (0.99 - 1)
    tau_xs = 0.5 * (alpha_xs + beta_xs)**-1
    xs_infinity = (1 + np.exp(V - 19.9 - 12.7))**-0.5
    dxs_dt = (xs_infinity - xs) / tau_xs

    # Component: L-type Ca channel (i_Ca_L)
    i_Ca_L = Cm * g_Ca_L * f * f_Ca * (V - 65)

    # Component: L-type Ca channel d gate
    d_infinity = (1 + np.exp(V + 10 - 8))**-1
    tau_d = np.where(abs(V + 10) < 1e-10, 4.5791 + np.exp(V + 10 - 6.24),
                   1 - np.exp(V + 10 - 6.24) * (0.035 * (V + 10) / (1 + np.exp(V + 10 - 6.24)))
    dd_dt = (d_infinity - d) / tau_d

    # Component: L-type Ca channel f gate
    f_infinity = np.exp(-(V + 28) / 6.91) + np.exp(-(V + 28) / 6.9)
    tau_f = 9 * ((0.0197 * np.exp(-(0.03372) * (V + 10)**2) + 0.02)**-1
    df_dt = (f_infinity - f) / tau_f

    # Component: L-type Ca channel f_Ca gate
    f_Ca_infinity = (1 + Ca_i / 0.00035)**-1
    tau_f_Ca = 2
    df_Ca_dt = (f_Ca_infinity - f_Ca) / tau_f_Ca

    # Component: sodium-potassium pump (i_NaK)
    f_NaK = (1 + 0.1245 * np.exp(-0.1 * F * V * R * T) + 0.0365 * sigma * np.exp(-F * V * R * T))**-1
    i_NaK = Cm * i_NaK_max * f_NaK / (1 + (Km_Na_i / Na_i)**1.5 * (K_o / K_o))

    # Component: background currents
    i_B_Na = Cm * g_B_Na * (V - E_Na)
    i_B_Ca = Cm * g_B_Ca * (V - E_Ca)
    i_B_K = Cm * g_B_K * (V - E_K)

    # Component: Na-Ca exchanger current (i_NaCa)
    i_NaCa = Cm * I_NaCa_max * (np.exp(gamma * F * V * R * T) * Na_i*3 * Ca_o - np.exp((gamma - 1) * F * V * R * T) * Na_o3 * Ca_i) / ((K_mNa3 + Na_o*3) * (K_mCa + Ca_o) * (1 + K_sat * np.exp((gamma - 1) * F * V * R * T)))

    # Component: sarcolemmal calcium pump current (i_CaP)
    i_CaP = Cm * i_CaP_max * Ca_i / (0.0005 + Ca_i)

    # Component: Ca release current from JSR (i_rel)
    i_rel = K_rel * u**2 * v * w * (Ca_rel - Ca_i)

    # Component: Ca release current from JSR u gate
    tau_u = 8
    u_infinity = (1 + np.exp(-(Fn - 3.4175e-13) / (13.67e-16)))**-1
    du_dt = (u_infinity - u) / tau_u

    # Component: Ca release current from JSR v gate
    tau_v = 1.91 + 2.09 * (1 + np.exp(-(Fn - 3.4175e-13) / (13.67e-16)))**-1
    v_infinity = 1 - ((1 + np.exp(-(Fn - 6.835e-14) / (13.67e-16)))**-1)
    dv_dt = (v_infinity - v) / tau_v

    # Component: Ca release current from JSR w gate
    tau_w = np.where(abs(V - 7.9) < 1e-10, 60.21 / 1.3,
                   6 / (1 - np.exp(-(V - 7.9) / 5)) * (1 + 0.3 * np.exp(-(V - 7.9) / 5) / (V - 7.9))
    w_infinity = 1 - ((1 + np.exp(-(V - 40) / 17))**-1)
    dw_dt = (w_infinity - w) / tau_w

    # Component: transfer current from NSR to JSR (i_tr)
    i_tr = (Ca_up - Ca_rel) / tau_tr

    # Component: Ca uptake current by the NSR (i_up)
    i_up = I_up_max / (1 + K_up / Ca_i)

    # Component: Ca leak current by the NSR (i_up_leak)
    #Ca_leak_current_by_the_NSR
    i_up_leak = I_up_max * Ca_up / (Ca_up_max + Ca_up)

       # Ca_buffers
    B1 = 2 * i_NaCa - (i_CaP + i_Ca_L + i_B_Ca) * 2 * Vi * F + Vup * (i_up_leak - i_up) + i_rel * Vrel * Vi
    B2 = 1 + Trpn_max * Km_Trpn * (Ca_i + Km_Trpn) ** 2 + Cmdn_max * Km_Cmdn * (Ca_i + Km_Cmdn) ** 2

    # Intracellular_ion_concentrations
    dCa_up = i_up - (i_up_leak + i_tr * Vrel * Vup)
    dCa_rel = (i_tr - i_rel) / (1 + Csqn_max * Km_Csqn * (Ca_rel + Km_Csqn) ** 2) - i_up

    return np.array([dV, dm, dh, dj, doa, doi, dua, dui, dxr, dxs, dd, df, df_Ca, dsigma, dF, dNa_i, dNa1o, dK_i, dK1o, dCa_i, di_NaCa_max, dgamma, dKm_Nai, dKm_Ca21o, di_Na, dE_Na, dalpha_m, dbeta_m, dm_inf, dtau_m, dalpha_h, dbeta_h, dh_inf, dtau_h, dalpha_j, dbeta_j, dj_inf, dtau_j, dE_K, di_K1, dalpha_oa, dbeta_oa, dtau_oa, doa_infinity, dalpha_oi, dbeta_oi, dtau_oi, doi_infinity, dg_Kur, di_Kur, dalpha_ua, dbeta_ua, dtau_ua, dua_infinity, dalpha_ui, dbeta_ui, dtau_ui, dui_infinity, dalpha_xr, dbeta_xr, dtau_xr, dxr_infinity, dalpha_xs, dbeta_xs, dtau_xs, dxs_infinity, di_Ca_L, dd_infinity, dtau_d, df_infinity, dtau_f, df_Ca_infinity, dtau_f_Ca, dE_Ca, di_B_Na, di_B_Ca, di_B_K, di_NaK, di_up, di_up_leak, dCa_CMDN, dKm_CMDN, dCa_TRPN, dKm_TRPN, dCa_CSQN, dKm_CSQN, dVi, dVrel, dVup])

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

# Define the initial conditions and time array
w0 = [-87.0,0.0,0.75,0.75,0.0,0.0,0.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.00007,0.5,0.5,18.0,145.0]
t = np.linspace(0, end_time, num_points)  # Adjust end_time and num_points accordingly

# Call the runge_kutta_4 function to solve the ODEs
results = runge_kutta_4(model, w0, t)
plt.figure(figsize=(12, 8))
plt.plot(t, results[:,0], linewidth=2 ,label= 'Euler')
plt.show()