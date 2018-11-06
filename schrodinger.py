#schrodinger.py
#author: Daniel Martin
#date created: 10/10/2018

"""
This python file solves Schrodinger's equation using
the shooting method.
"""
__author__ = 'Daniel Martin'
__version__ = '0.1'

"""
Imports!! Work out how to do
try: 
	blah
except ImportError:
	raise Error("message")

for importing numpy (if necessary) and any other imports. 
"""
try:
	import numpy as np
	import matplotlib.pyplot as plt
except ImportError:
	raise Error('The schrodinger class requires numpy and matplotlib to be installed.\
		Please install these modules before continuing.')



class schrodinger: #http://jakevdp.github.io/blog/2012/09/05/quantum-python/
	"""
	Class to implement numerical solution of TDSE 
	using the shooting method for an arbitrary potential.
	"""
	
	def __init__(self, x, V, hbar=1, m=1, e0=1):
		"""
		Parameters

		x : array of floats
			array of length N of evenly spaced spatial coordinates
		V : array of floats
			array of length N defining potential for x
		hbar : float
			value of the reduced Planck's constant (default = 1) 
		m : float
			value of particle's mass (default = 1)
		e0 : float
			value of electronic charge (default = 1)
		"""

		# Set initial conditions
		self.hbar = hbar
		self.m = m
		self.e0 = e0
		self.x = x
		self.V = V
		self.N = self.x.size
		self.dx = abs(self.x[1] - self.x[0])

		# Assert inputs of correct type and size
		assert isinstance(self.x, np.ndarray), 'x must be a numpy array'
		assert isinstance(self.V, np.ndarray), 'V must be a numpy array'
		assert self.x.shape == (self.N,), 'x and V must have same size'
		assert self.V.shape == (self.N,), 'x and V must have same size'
		assert isinstance(self.hbar, (int, float))
		assert isinstance(self.m, (int, float))
		assert isinstance(self.e0, (int, float))

	""" 		METHOD AS OUTLINED IN HARRISON, COMPUTATIONAL METHODS IN PHYSICS, CHEMISTRY AND BIOLOGY
	def psi_at_inf(self, E):
		psi0 = 0
		psi1 = 1

		for i in range(self.N):
			psi2 = (2 * self.m * (self.dx / self.hbar)**2 * \
				   (self.e0 * (self.x[i] / 100e-10)**2 - E) + 2) * psi1 - psi0
			psi0 = psi1
			psi1 = psi2

		return psi2

	def shooting(self):
		dE = 1e-3 * self.e0
		eigvals = []

		Y1 = self.psi_at_inf(0)
		E = 0

		for i in range(self.N):
			E += dE
			Y2 = self.psi_at_inf(E)

			if Y1*Y2 < 0:
				eigvals.append((abs(Y1)*dE/(abs(Y1)+abs(Y2))+E-dE)/(1e-3*e0))

		return eigvals
	"""

	


	"""
	def plot_eigenfunctions(self, samefig = True):
		assert isinstance(samefig, bool), "samefig must be a boolean"

		# Test if matplotlib.pyplot is imported as plt
		try:
			plt.figure("test")
		except NameError:
			try:
				import matplotlib.pyplot as plt
			except ImportError:
				print("In order to plot the eigenfunctions, matplotlib must be \
						installed on the system")
		plt.close("test")

		# Test to see if figure already exists
		if samefig == True and plt.fignum_exists(1) == True:
			plt.figure(1)
		else:
			plt.figure()

		# Plot potential
		plt.plot(self.x, self.V, 'k', label="potential")

		# Plot wavefunctions										******* DO THIS ******
		


	def plot_eigenenergies(self, samefig = True):
		assert isinstance(samefig, bool), "samefig must be a boolean"

		# Test if matplotlib.pyplot is imported as plt
		try:
			plt.figure("test")
		except NameError:
			try:
				import matplotlib.pyplot as plt
			except ImportError:
				print("In order to plot the eigenenergies, matplotlib must be \
						installed on the system")
		plt.close("test")

		# Test to see if figure already exists
		if samefig == True and plt.fignum_exists(1) == True:
			plt.figure(1)
		else:
			plt.figure()

		# Plot potential
		plt.plot(self.x, self.V, 'k', label="potential")

		# Plot wavefunctions										******* DO THIS ******
		


	def plot(self):
		plot_eigenenergies()
		plot_eigenfunctions()
		plt.show()
	"""


x = np.linspace(0, 6e-9, 1000)
v = np.zeros(1000)
a = schrodinger(x, v)

eig = a.shooting()

print(a.N)