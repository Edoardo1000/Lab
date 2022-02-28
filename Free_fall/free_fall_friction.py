import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def model(t, g, a, v0, y0):
	return y0 - g*t / a + (v0/a + g/a**2)*(1-np.exp(-a*t))  

#################### Medatada ####################  
filename = "data2/biliardo.txt"
fileTitle = "Palla da biliardo" + " (con attrito)"
fileO = "grafici/billiard_refined.png"
#################### Data loading ####################  

data = np.loadtxt(filename)
y = data[1:,1]
n_frames = len(y) 
sigma_y = np.linspace(data[0,0]/100, data[0,1]/100 ,n_frames)
t = data[1:,0]


#################### Fitting ####################  
popt, pcov = curve_fit(model, t, y, sigma = sigma_y, maxfev=10000, bounds=((0,0,-np.inf,-np.inf),(20,0.1,10,20)))
g, a, v0, y0 = popt

sigma_g = np.sqrt(pcov.diagonal()[0])
sigma_a = np.sqrt(pcov.diagonal()[1])
print(f"g = {g} | sigma_g = {sigma_g} | a = {a} | sigma_a = {sigma_a} | v0 = {v0} | y0 = {y0}")

#################### Graphing ####################  
fig, axs = plt.subplots(2)

axs[0].errorbar(t, y, sigma_y, fmt='.')
axs[0].grid(visible=True, which="both", ls="dashed", color="gray")
axs[0].set_ylabel("altezza [m]")
axs[0].plot(t, model(t, g, a, v0, y0))
axs[0].set_title(fileTitle)
axs[1].errorbar(t, y - model(t, g, a, v0, y0), sigma_y, fmt=".")
axs[1].set_ylabel("residui [m]")
axs[1].set_xlabel("tempo [s]")
axs[1].grid(visible=True, which="both", ls="dashed", color="gray")
axs[1].plot(t, np.zeros(len(t)))

plt.savefig(fileO, dpi=400)

#################### Chi squared ####################  
zeta = np.square(np.divide(y - model(t, g, a, v0, y0), sigma_y))/(n_frames-4)
chi_normalized = np.sum(zeta)
print(f"chi_normalized = {chi_normalized}")
