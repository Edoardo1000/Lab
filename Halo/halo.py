import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('point-coordinates.txt', unpack=True, delimiter=',')
x_s = x[0]
y_s = y[0]

x = x[1:]
y = y[1:]
n = len(x)
u = x - x.mean()
v = y - y.mean()

s_u = np.sum(u)
s_uu = (u**2).sum()
s_uuu = (u**3).sum()

s_v = np.sum(v)
s_vv = (v**2).sum()
s_vvv = (v**3).sum()

s_uv = (u*v).sum()
s_uuv = (u*u*v).sum()
s_uvv = (u*v*v).sum()
D = 2*(s_uu*s_vv-s_uv**2)

u_c = (s_vv*(s_uuu+s_uvv)-s_uv*(s_vvv+s_uuv))/D
v_c = (s_uu*(s_vvv+s_uuv)-s_uv*(s_uuu+s_uvv))/D
r = np.sqrt(u_c**2 + v_c**2 + (s_uu + s_vv)/n)

x_c = u_c + x.mean()
y_c = v_c + y.mean()

print(x_c,y_c,r)
######################### Plotting #########################

plt.figure("Fit circolare")
x_space = np.linspace(300,850)
y_space = np.linspace(175,600)

X, Y = np.meshgrid(x_space,y_space)
F = (X-x_c)**2 + (Y-y_c)**2 - r**2
plt.contour(X, Y, F, [0])

plt.scatter(x,y,label="punti misurati")
plt.scatter(x_s,y_s,label='posizione del centro')
plt.title("Fit dell'alone lunare")
plt.legend()

plt.grid(ls='dashed',color='grey',which='both')

plt.savefig('halo.jpg',dpi=500)