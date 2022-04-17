import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def model(t, omega_0, omega_k, gamma, A, B, phi_1, phi_2, C):
	return C + np.exp(-gamma * t)*(A * np.cos( np.sqrt(omega_0**2 - gamma**2)*t + phi_1) + B * np.cos( np.sqrt(omega_0**2 + 2 * omega_k ** 2 - gamma**2)*t + phi_2 )  )

t_a, x_a, t_b, x_b = np.loadtxt("./dati/combination_1.txt", unpack=True)

x_a=x_a/100
x_b=x_b/100

sigma_a = 0.005/np.sqrt(12) * x_a
sigma_b = 0.005/np.sqrt(12) * x_b

#################### Fitting ####################
popt_a, pcov_a = curve_fit(model, t_a, x_a, sigma=sigma_a, p0=(4, 2, 0.1, 10, 10, 0, 0, 200), maxfev=10000, absolute_sigma=True)

omega_0_a, omega_k_a, gamma_a, A_a, B_a, phi_1_a, phi_2_a, C_a = popt_a
sigma_omega_0_a, sigma_omega_k_a, sigma_gamma_a, sigma_A_a, sigma_B_a, sigma_phi_1_a, sigma_phi_2_a, sigma_C_a = np.sqrt(pcov_a.diagonal())


popt_b, pcov_b = curve_fit(model, t_b, x_b, sigma=sigma_b, p0=(4, 2, 0.1, 200, 200, 0, 0, 200), maxfev=10000, absolute_sigma=True)

omega_0_b, omega_k_b, gamma_b, A_b, B_b, phi_1_b, phi_2_b, C_b = popt_b
sigma_omega_0_b, sigma_omega_k_b, sigma_gamma_b, sigma_A_b, sigma_B_b, sigma_phi_1_b, sigma_phi_2_b, sigma_C_b = np.sqrt(pcov_b.diagonal())

print("omega_a = %.5f+/-%.5f" % (omega_0_a, sigma_omega_0_a),"omega_b = %.5f+/-%.5f" % (omega_0_b, sigma_omega_0_b))
print("omegak_a = %.5f+/-%.5f" % (omega_k_a, sigma_omega_k_a),"omegak_b = %.5f+/-%.5f" % (omega_k_b, sigma_omega_k_b))
print("gamma_a = %.5f+/-%.5f" % (gamma_a, sigma_gamma_a),"gamma_b = %.5f+/-%.5f" % (gamma_b, sigma_gamma_b))
#################### Plotting ####################

plt.subplot(2,1,1)
plt.errorbar(t_a, x_a, sigma_a, fmt='.')
plt.errorbar(t_b, x_b, sigma_b, fmt='.')
plt.plot(t_a, model(t_a, omega_0_a, omega_k_a, gamma_a, A_a, B_a, phi_1_a, phi_2_a, C_a))
plt.plot(t_b, model(t_b, omega_0_b, omega_k_b, gamma_b, A_b, B_b, phi_1_b, phi_2_b, C_b))
plt.title("un oscillatore inclinato")
plt.ylabel('ampiezza [ua]')
plt.grid(ls='dashed',which='both')
plt.xlim([0,45])

plt.subplot(2,1,2)
plt.errorbar(t_a, x_a - model(t_a, omega_0_a, omega_k_a, gamma_a, A_a, B_a, phi_1_a, phi_2_a, C_a),sigma_a,fmt='.')
plt.errorbar(t_b, x_b - model(t_b, omega_0_b, omega_k_b, gamma_b, A_b, B_b, phi_1_b, phi_2_b, C_b),sigma_b,fmt='.')
plt.grid(ls='dashed',which='both')
plt.xlabel('tempo [s]')
plt.ylabel('residui [ua]')
plt.xlim([0,45])

#################### Plotting ####################
plt.savefig('combination.jpg',dpi=500)
plt.show()    