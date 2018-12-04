# Analytic energy values for parabolic and linear well

import numpy as np
mev=1.60217656e-22; hbar=1.0545718e-34; q=1.60217656e-19; kb=1.38e-23; T=300.; eps0=8.85e-12; m_e=9.10938356e-31; Ang=1.0e-10

def linearE(n):
	global hbar, m_e, mev
	left = (hbar**2 / (2*0.014*m_e))**(1./3)
	c = 150*mev / 80e-9
	right = (1.5*np.pi*c*(n-0.25))**(2./3)
	return left*right / mev


def ParabolicE(n):
	global hbar, m_e, mev
	alpha = 2* 150*mev / (40e-9)**2
	omega = (alpha / (0.014*m_e))**0.5
	return  (n+0.5)*hbar*omega / mev

for i in range(1, 5):
	print(f'Parabola: E({i}) = {ParabolicE(i)}')

print('\n')

for i in range(1, 5):
	print(f'Linear: E({i}) = {linearE(i)}')