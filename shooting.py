import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import newton
plt.close("all")
"""
# FIRST ATTEMPT AT SHOOTING FROM WELL WIRES AND DOTS. GIVES EIGENFUNCTIONS WHEN EIGENENERGY PUT IN.
def shooting(x, V, E, m=1, hbar=1):
	dx = abs(x[1] - x[0])

	# decay = np.sqrt(2*m*(V - E)) / hbar
	psi = np.zeros(len(x))
	a = c = hbar**2 / (2*m*dx**2)
	b = hbar**2 / (m*dx**2) + V
	psi[0] = 0
	psi[1] = 1
	for i in range(2, len(x)):
		psi[i] = (1/c) * (E - b[i])*psi[i-1] - psi[i-2]
	return psi


x = np.linspace(0, 1, 1001)
V = np.zeros(len(x))

plt.figure()
for i in range(1, 5):
	psi = shooting(x, V, i**2 * np.pi**2/2) 
	plt.plot(x[::2], -psi[::2])
	print(max(psi))

# plt.figure()
# plt.plot(x, psi) # Gives block
# plt.plot(x[::2], -psi[::2]) # Gives smooth line

"""

"""
def RK4(SE, psi0, x):
	N = len(x)
	dx = abs(x1[1] - x[0])
	psi = np.ones(N) * psi0    #Sets initial value of psi to psi0
	for i in range(N-1):
		k1 = dx * SE(x[i], psi[i])
		k2 = dx * SE(x[i] + dx/2.0, psi[i] + k1/2.0)
		k3 = dx * SE(x[i] + dx/2.0, psi[i] + k2/2.0)
		k4 = dx * SE(x[i] + dx, psi[i] + k3)
		psi[i+1] = psi[i] + (k1 + 2*k2 + 2*k3 + k4) / 6.0
	return psi
"""



# Harrison Computational Methods for Physics, Bio, Chem
def psi_at_inf(E):
	hbar = 1.05459e-34
	m = 9.109534e-31
	e0 = 1.602189e-19
	dx = 1e-10
	psi0 = 0
	psi1 = 1
	x = dx
	while (x < 100e-10):
		psi2 = (2 * m * (dx / hbar) * (dx / hbar) * (e0 * (x / 100e-10) * (x / 100e-10) - E) + 2) * psi1 - psi0
		psi0 = psi1
		psi1 = psi2
		x += dx
	return psi2

def shooting():
	hbar = 1.05459e-34
	m = 9.109534e-31
	e0 = 1.602189e-19
	dE = 1e-3 * e0
	eigenergy = []

	Y1 = psi_at_inf(0)
	E = dE
	while (E < e0):
		Y2 = psi_at_inf(E)
		if (Y1*Y2 < 0):
			eig = (abs(Y1) * dE / (abs(Y1) + abs(Y2)) + E - dE) / (1e-3 * e0)
			eigenergy.append(eig)
		Y1 = Y2
		E += dE
	return eigenergy




x = np.linspace(0, 1, 1000)
V = np.ones(len(x))
V[450:550] = 0
eig = shooting()
for vals in eig:
	print(f'energy eigval at {vals:.2f} meV')
# print(f'E = {E}')

# plt.figure()
# plt.plot(x, psi)
# plt.show()