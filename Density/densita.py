import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# models for the fit
def line(x, m, q):
    return m * x + q

def power_law(x, norm, index):
    return norm * (x**index)

###############################  Brass  ###############################
# experimental data
m_brass = np.array([10.490, 24.478, 16.522, 39.954])
sigma_m_brass = np.full(m_brass.shape, 0.001)

V_brass = np.array([1247.8, 2906.1, 2004.6, 4136.7])
sigma_V_brass = np.array([2.6, 7.0, 4.1, 7.7])

# fit of the data, lambda is the inverse of the density
popt_brass, pcov_brass = curve_fit(line, m_brass, V_brass, sigma=sigma_V_brass, absolute_sigma=True)
lambda_hat_brass, q_hat_brass = popt_brass
sigma_lambda_brass, sigma_q_brass = np.sqrt(pcov_brass.diagonal())
residue_brass = V_brass - line(m_brass, lambda_hat_brass, q_hat_brass)

rho_hat_brass = 10**6/lambda_hat_brass
sigma_rho_brass = sigma_lambda_brass/lambda_hat_brass**2 * 10**6




# data plot
plt.figure("mass-volume graph (brass)")
plt.subplot(2,1,1)
plt.errorbar(m_brass, V_brass, sigma_V_brass, sigma_m_brass, fmt='.')
plt.title("Volume dell'ottone in funzione della massa")
plt.grid(which="both")
plt.ylabel("Volume $[\mathrm{mm}^3]$")

# fit and residue plot
x = np.linspace(0,42,200)
plt.plot(x, line(x,lambda_hat_brass, q_hat_brass))

plt.subplot(2,1,2)
plt.errorbar(m_brass, residue_brass, sigma_V_brass, sigma_lambda_brass, fmt='.')
plt.plot(x,np.zeros(200))
plt.grid(which="both")
plt.xlabel("Mass [g]")
plt.ylabel("Volume $[\mathrm{mm}^3]$")

plt.savefig("brass.png")
############################## Alluminium ##############################
# experimental data
m_al = np.array([5.848, 1.452, 4.86, 7.770])
sigma_m_al = np.full(m_al.shape, 0.001)

V_al = np.array([2162.4, 543.18, 1813.8, 2899.5])
sigma_V_al = np.array([3.8, 1.8, 2.7, 3.6])

# fit of the data, lambda is the inverse of the density
popt_al, pcov_al = curve_fit(line, m_al, V_al, sigma=sigma_V_al, absolute_sigma=True)
lambda_hat_al, q_hat_al = popt_al
sigma_lambda_al, sigma_q_al = np.sqrt(pcov_al.diagonal())
residue_al = V_al - line(m_al, lambda_hat_al, q_hat_al)

rho_hat_al = 10**6/lambda_hat_al
sigma_rho_al = sigma_lambda_al/lambda_hat_al**2 * 10**6

# data plot
plt.figure("mass-volume graph (alluminium)")
plt.subplot(2,1,1)
plt.errorbar(m_al, V_al, sigma_V_al, sigma_m_al, fmt=".")
plt.title("Volume dell'alluminio in funzione della massa")
plt.ylabel("Volume $[\mathrm{mm}^3]$")
plt.grid(which="both")

# fit and residue plot
x = np.linspace(0,8,200)
plt.plot(x, line(x, lambda_hat_al,q_hat_al))

plt.subplot(2,1,2)
plt.errorbar(m_al, residue_al, sigma_V_al, sigma_m_al, fmt='.')
plt.plot(x,np.zeros(200))
plt.grid(which="both")
plt.xlabel("Mass [g]")
plt.ylabel("Volume $[\mathrm{mm}^3]$")

plt.savefig("alluminium.png")

print(f"rho_brass = {rho_hat_brass} | q_brass = {q_hat_brass} | sigma_rho_brass = {sigma_rho_brass} | sigma_q_brass = {sigma_q_brass}\n")
print(f"rho_al = {rho_hat_al} | q_al = {q_hat_al} | sigma_rho_al = {sigma_rho_al} | sigma_q_al = {sigma_q_al}\n")



plt.show()