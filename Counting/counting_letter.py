import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, norm, binom


#################### Collecting data ####################

lenghts = []
line_lenghts = []
l = []

with open('commedia.txt') as testo:
	lines = testo.readlines()
	lines_better = []

	for line in lines:
		if len(line) > 3:
			line_fixed = line.replace("\n","").lstrip().lower().replace('à','a').replace('è','e')
			lines_better.append(line_fixed)
			lenghts.append(len(line_fixed))
			l.append(line_fixed.count("a"))
			line_lenghts.append(len(line_fixed))

l = np.array(l)
line_lenghts = np.array(line_lenghts)
print(l)
print(len(line_lenghts))

binning = np.arange(l.min()-0.5,l.max()+1.5)
o,_,_ = plt.hist(l, bins=binning, rwidth=0.25, label="Conteggi")
p = (l.sum()+0.0)/line_lenghts.sum()
N=line_lenghts.mean()
print(p,N)
#################### Plotting ####################

k=np.arange(l.min(),l.max()+1)

e_binom = binom.pmf(k,N,p)*len(line_lenghts)
e_poisson = poisson.pmf(k,p*N)*len(line_lenghts)


chi_binom = ((o-e_binom)**2/e_binom).sum()
chi_poisson = ((o-e_poisson)**2/e_poisson).sum()
dof_binom = len(k) - 2 - 1
dof_poisson = len(k) - 1 - 1

print(f"Chi quadro per la binomiale:{chi_binom} \ {dof_binom}")
print(f"Chi quadro per la poissoniana:{chi_poisson} \ {dof_poisson}")

plt.bar(k-0.3,e_binom,width=0.25, color="#ff7f0e",label="Binomiale")
plt.bar(k+0.3,e_poisson,width=0.25, color="#2ca02c",label="Poisson")


plt.legend()
plt.grid(ls="dashed",color="grey",which="both")

plt.show()
