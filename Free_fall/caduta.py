import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parabola(x, a, b, c):
    return a*x**2 + b*x + c


#################### Medatada ####################  
filename = "data/billiard.txt"
fileTitle = "Billiard Ball"
fileO = "graphs/billiard.png"
#################### Data loading ####################  

data = np.loadtxt(filename)
y = data[2:,1]
n_frames = int(data[1,1])


time_elapsed = (n_frames)/240.0
sigma_y = np.linspace(data[0,0]/100, data[0,1]/100 ,n_frames+1)
t = np.linspace(0,time_elapsed, n_frames+1)


#################### Fitting ####################  
popt, pcov = curve_fit(parabola, t, y, sigma = sigma_y)
g, v0, h = popt

sigma_g = np.sqrt(pcov.diagonal()[0])
print(f"g = {2*g} | sigma_g = {2*sigma_g}")

#################### Graphing ####################  
fig, axs = plt.subplots(2)

axs[0].errorbar(t, y, sigma_y, fmt='.')
axs[0].grid(visible=True, which="both", ls="dashed", color="gray")
axs[0].set_ylabel("height [m]")
axs[0].plot(t, parabola(t, g, v0, h))
axs[0].set_title(fileTitle)
axs[1].errorbar(t, y - parabola(t, g, v0, h), sigma_y, fmt=".")
axs[1].set_ylabel("residues [m]")
axs[1].grid(visible=True, which="both", ls="dashed", color="gray")

plt.savefig(fileO)

#################### Chi squared ####################  
zeta = np.square(np.divide(y - parabola(t, g, v0, h), sigma_y))/(n_frames-3)
chi_normalized = np.sum(zeta)
print(f"chi_normalized = {chi_normalized}")