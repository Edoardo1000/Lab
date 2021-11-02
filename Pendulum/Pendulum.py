import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Calculates the Standard Deviation of the Mean of a list x
def SDM(x):
  m = np.mean(x)
  
  sum = 0
  for e in x:
    sum += (e-m)**2
  sum /= (len(x)*(len(x)-1))
  sum = np.sqrt(sum)

  return sum

# The function used for the fit, l is the free parameter
def period_model(d, l, g=9.81):
  return 2.0 * np.pi * np.sqrt((l**2 / 12.0 + d**2) / (g*d))


# Opening, reading and closing the csv file with data
csv_file = open("Pendolo_Fisico_Dati.csv", mode="r")

reader = csv.reader(csv_file, delimiter=',')
header = next(reader)
data = []
for row in reader:
  data.append(row)

csv_file.close()

# Putting all the data in numpy arrays
data_array = np.array(data)
T = []
sigma_T = []
for i in range(data_array.shape[0]):
  row = data_array[i,1:].astype(float)
  mean = np.mean(row)
  T.append(mean)
  sigma_T.append(SDM(row))

d = data_array[:,0].astype(float)
sigma_d = np.full(d.shape, 0.0011)
T = np.array(T).astype(float)
sigma_T = np.array(sigma_T)


# fitting
popt, pcov = curve_fit(period_model, d, T, sigma=sigma_T)
l_hat = popt[0]
sigma_l = np.sqrt(pcov[0, 0])
print("l =",l_hat)
print("sigma l =", sigma_l)


# Graph of model and residues
fig, axs = plt.subplots(2)
axs[0].errorbar(d, T, sigma_T, sigma_d, fmt='.')
x = np.linspace(0.025, 0.5, 100)
axs[0].plot(x, period_model(x, l_hat))
axs[0].grid(which='both', ls='dashed', color='gray')

Sigma_l = np.full(d.shape, sigma_l)
axs[1].errorbar(d, T - period_model(d, l_hat), sigma_T, fmt='.')
axs[1].grid(which='both', ls='dashed', color='gray')
axs[1].plot(x, np.full(x.shape, 0))

axs[1].set_xlabel("Distanza del centro di massa [m]")
axs[0].set_ylabel("Periodo [s]")
axs[1].set_ylabel("Residuo [s]")

plt.savefig("Pendolo.jpg")
plt.show()
