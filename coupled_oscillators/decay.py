import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




def model(y, x_0, tau, C, B):
	
	return x_0 * np.exp(-y /tau) + C + B*np.exp(-3*y/tau)

	

t, x, t_b, x_b = np.loadtxt("dati/damped_1.txt", unpack=True)

t_max = []
x_max = []

for i in range(1, len(t)-1):
	if x[i] >= x[i-1] and x[i] >= x[i+1] and x[i]>450:
		t_max.append(t[i])
		x_max.append(x[i])

t_max = np.array(t_max)
x_max = np.array(x_max)

x_max = x_max/100
sigma_x = x_max * 0.001

popt, pcov = curve_fit(model, t_max, x_max, p0=(x_max[0], 40, 443, 0), bounds=((0, 0, -1000, -500),(10000, 80, 600, 500)))
x_0, tau, C, B= popt
sigma_x_0, sigma_tau, sigma_C, sigma_B = np.sqrt(pcov.diagonal())
print("tau = %.5f +/- %.5f" % (tau, sigma_tau))

#################### Plotting #################### 
plt.subplot(2,1,1)

t_space = np.linspace(0, t[-1]+5, 100)
plt.errorbar(t_max, x_max, sigma_x, fmt='.')
plt.plot(t_space, model(t_space, x_0, tau, C, B))
plt.grid(ls='dashed', color='grey', which='both')
plt.ylabel('ampiezza [ua]')
plt.title("decadimento con correzione")

plt.subplot(2,1,2)
plt.errorbar(t_max, x_max - model(t_max, x_0, tau, C, B), sigma_x, fmt='.')
plt.grid(ls='dashed', color='grey', which='both')
plt.plot(t_space, np.zeros(len(t_space)))
plt.savefig("figure.jpg", dpi=500)
plt.ylabel('residui [ua]')
plt.xlabel('tempo [s]')
plt.show()