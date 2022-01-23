from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# models for the fit
def line(x, m, q):
    return m * x + q

def power_law(x, norm, index):
    return norm * (x**index)

############################## Allumium Direct ##############################
# experimental data
m_al = np.array([5.848, 1.452, 4.86, 7.770])
sigma_m_al = np.full(m_al.shape, 0.001)

V_al = np.array([2162.4, 543.18, 1813.8, 2899.5])
sigma_V_al = np.array([3.8, 1.8, 2.7, 3.6])

popt_al_direct, pcov_al_direct = curve_fit(line, V_al, m_al, sigma=sigma_m_al, absolute_sigma=True)
rho_al, q_al = popt_al_direct
sigma_rho_al, sigma_q_al = np.sqrt(pcov_al_direct.diagonal())
residue_al = m_al - line(V_al, rho_al, q_al)

# data plot
plt.figure("direct mass-volume graph (al)")
plt.subplot(2,1,1)
plt.errorbar(V_al, m_al, sigma_m_al, sigma_V_al, fmt='.')
plt.title("Massa dell'alluminio in funzione del volume")
plt.grid(which="both")
plt.ylabel("mass [g]")

# fit and residue plot
x = np.linspace(0,4200,200)
plt.plot(x, line(x,rho_al, q_al))

plt.subplot(2,1,2)
plt.errorbar(V_al, residue_al, sigma_m_al, sigma_rho_al, fmt='.')
plt.plot(x,np.zeros(200))
plt.grid(which="both")
plt.xlabel("volume $[\mathrm{mm}^3]$")
plt.ylabel("Mass [g]")

plt.savefig("alluminium_direct.png")
print(f"rho_al = {rho_al*10**6} | q_al = {q_al*1000} | sigma_rho_al = {sigma_rho_al*10**6} | sigma_q_al = {sigma_q_al}\n")



# experimental data
m_brass = np.array([10.490, 24.478, 16.522, 39.954])
sigma_m_brass = np.full(m_brass.shape, 0.001)

V_brass = np.array([1247.8, 2906.1, 2004.6, 4136.7])
sigma_V_brass = np.array([2.6, 7.0, 4.1, 7.7])

############################## Brass Direct ##############################
popt_brass, pcov_brass = curve_fit(line, V_brass, m_brass, sigma=sigma_m_brass, absolute_sigma=True)
rho_brass, q_brass = popt_brass
sigma_rho_brass, sigma_q_brass = np.sqrt(pcov_brass.diagonal())
residue_brass = m_brass - line(V_brass, rho_brass, q_brass)

# data plot
plt.figure("direct mass-volume graph (brass)")
plt.subplot(2,1,1)
plt.errorbar(V_brass, m_brass, sigma_m_brass, sigma_V_brass, fmt='.')
plt.title("Massa dell'ottone in funzione del volume")
plt.grid(which="both")
plt.ylabel("mass [g]")

# fit and residue plot
x = np.linspace(0,4200,200)
plt.plot(x, line(x,rho_brass, q_brass))

plt.subplot(2,1,2)
plt.errorbar(V_brass, residue_brass, sigma_m_brass, sigma_rho_brass, fmt='.')
plt.plot(x,np.zeros(200))
plt.grid(which="both")
plt.xlabel("volume $[\mathrm{mm}^3]$")
plt.ylabel("Mass [g]")
plt.savefig("brass_direct.png")

print(f"rho_brass = {rho_brass*10**6} | q_brass = {q_brass*1000} | sigma_rho_brass = {sigma_rho_brass*10**6} | sigma_q_brass = {sigma_q_brass}\n")
