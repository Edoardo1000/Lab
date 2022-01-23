import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# models for the fit
def line(x, m, q):
    return m * x + q

def power_law(x, norm, index):
    return norm * (x**index)


# Values for the spheres
# mass
m = np.array([8.359, 11.894, 18.913, 24.768, 44.837])
sigma_m = np.full(m.shape, 0.001)
# diameter
d = np.array([12.70, 14.29, 16.68, 18.25, 22.23])
sigma_d = np.full(d.shape, 0.01)
# volume
V = 4.0 / 24.0 * np.pi * d**3.0
sigma_V = V * 3.0 * sigma_d / d

popt_sphere, pcov_sphere = curve_fit(power_law, d, m, sigma=sigma_m, absolute_sigma=True)
norm_hat, index_hat = popt_sphere
sigma_norm, sigma_index = np.sqrt(pcov_sphere.diagonal())

x = np.linspace(5, 25, 100)
residue_sphere = m - power_law(d, norm_hat, index_hat)

plt.figure("mass-diameter graph (spheres)")
plt.subplot(2,1,1)
plt.errorbar(d, m, sigma_m, sigma_d, fmt='.')
plt.xlabel("diameter [mm]")
plt.ylabel("volume [$\mathrm{mm}^3$]")
plt.xscale('log')
plt.yscale('log')
plt.grid(which='both', ls='dashed')
plt.title("Massa delle sfere in funzione del diametro")
plt.plot(x, power_law(x, norm_hat, index_hat))

plt.subplot(2,1,2)
plt.errorbar(d, residue_sphere, sigma_m, fmt='.')
plt.plot(x,np.zeros(100))
plt.grid(which='both', ls='dashed')
plt.xlabel('diameter[mm]')
plt.ylabel('redisue [$\mathrm{mm}^3$]')

plt.savefig("spheres.png")

print(f"index = {index_hat} | norm = {norm_hat} | sigma_index = {sigma_index} | sigma_norm = {sigma_norm}")


############################## Mass-Volume ##############################
popt_line, pcov_line = curve_fit(line, m, V, sigma=sigma_V, absolute_sigma=True)
rho_hat_inverse, q_hat_inverse = popt_line
sigma_rho_hat_inverse, sigma_q_hat_inverse = np.sqrt(pcov_line.diagonal())
x = np.linspace(0, 45, 100)


rho_hat = 10**6/rho_hat_inverse
sigma_rho = sigma_rho_hat_inverse/rho_hat_inverse**2 * 10**6

plt.figure("volume-mass graph (spheres)")
plt.subplot(2,1,1)
plt.errorbar(m,V, sigma_V, sigma_m, fmt=".")
plt.title("Volume delle sfere in funzione della massa")
plt.ylabel('volume [$\mathrm{mm}^3]$')
plt.grid(which='both', ls='dashed')
plt.plot(x, line(x, rho_hat_inverse, q_hat_inverse))

residue = V - line(m, rho_hat_inverse, q_hat_inverse)

plt.subplot(2,1,2)
plt.errorbar(m, residue, sigma_V, fmt='.')
plt.xlabel('mass [g]')
plt.ylabel('volume [$\mathrm{mm}^3]$')
plt.plot(x, np.zeros(100))
plt.grid(which='both', ls='dashed')

print(f"inverse: rho_hat = {rho_hat} | sigma_rho = {sigma_rho} | q_hat_inverse = {q_hat_inverse} | sigma_q_hat_inverse = {sigma_q_hat_inverse}")
plt.savefig("steel.png")
############################## Volume-Mass ##############################

popt_line, pcov_line = curve_fit(line, V, m, sigma=sigma_m, absolute_sigma=True)
rho_hat, q_hat= popt_line
sigma_rho_hat, sigma_q_hat= np.sqrt(pcov_line.diagonal())
x = np.linspace(0, 6000, 100)


plt.figure("mass-volume graph (spheres)")
plt.subplot(2,1,1)
plt.title("Massa delle sfere in funzione del volume")
plt.errorbar(V,m, sigma_m, sigma_V, fmt=".")
plt.ylabel('mass [g]')
plt.grid(which='both', ls='dashed')
plt.plot(x, line(x, rho_hat, q_hat))

residue = m - line(V, rho_hat, q_hat)

plt.subplot(2,1,2)
plt.errorbar(V, residue, sigma_m, fmt='.')
plt.ylabel('mass [g]')
plt.xlabel('volume [$\mathrm{mm}^3]$')
plt.plot(x, np.zeros(100))
plt.grid(which='both', ls='dashed')

print(f"direct: rho_hat = {10**6*rho_hat} | sigma_rho = {sigma_rho} | q = {q_hat} | sigma_q = {sigma_q_hat}")
plt.savefig("steel_direct.png")

