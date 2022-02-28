import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parabola(x, a, b, c):
    return a*x**2 + b*x + c


#################### Medatada ####################  
filename = "data2/biliardo.txt"
fileTitle = "Palla da biliardo"
fileO = "grafici/billiard.png"
#################### Data loading ####################  

data = np.loadtxt(filename)
y = data[1:,1]
n_frames = len(y) 
sigma_y = np.linspace(data[0,0]/100, data[0,1]/100 ,n_frames)
t = data[1:,0]


#################### Fitting ####################  
popt, pcov = curve_fit(parabola, t, y, sigma = sigma_y)
g, v0, h = popt

sigma_g = np.sqrt(pcov.diagonal()[0])
print(f"g = {2*g} | sigma_g = {2*sigma_g}")

#################### Graphing ####################  
fig, axs = plt.subplots(2)

axs[0].errorbar(t, y, sigma_y, fmt='.')
axs[0].grid(visible=True, which="both", ls="dashed", color="gray")
axs[0].set_ylabel("altezza [m]")
axs[0].plot(t, parabola(t, g, v0, h))
axs[0].set_title(fileTitle)
axs[1].errorbar(t, y - parabola(t, g, v0, h), sigma_y, fmt=".")
axs[1].set_ylabel("residui [m]")
axs[1].set_xlabel("tempo [s]")
axs[1].plot(t, np.zeros(len(t)))
axs[1].grid(visible=True, which="both", ls="dashed", color="gray")

plt.savefig(fileO, dpi=400)

#################### Chi squared ####################  
zeta = np.square(np.divide(y - parabola(t, g, v0, h), sigma_y))/(n_frames-4)
chi_normalized = np.sum(zeta)
print(f"chi_normalized = {chi_normalized}")
