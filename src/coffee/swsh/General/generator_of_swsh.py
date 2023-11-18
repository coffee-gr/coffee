#!/usr/bin/env python
# encoding: utf-8
"""
Here we construct the spin weighted-spherical harmonics based on the formulas 
given in:
Goldberg, Joshua N., et al. Spin-s Spherical Harmonics and eths. 
Journal of Mathematical Physics 8.11 (1967): 2155-2161.
"""

from coffee.settings import be
import math
import sys
import sympy as sm

def Binomial(x,y):		
	a = math.factorial(x)
	if (y < 0):
		return 0
	elif (x - y) < 0:
		return 0
	else:
		b = math.factorial(y)
		c = math.factorial(x - y)
		return sm.Rational(a, (b * c))

def Formula_of_Yslm(s,l,m):
	 		    
	theta, phi = sm.symbols('theta phi')

	Numerical_part = sm.Rational(math.factorial(l+m) * math.factorial(l-m) \
				* (2*l+1), 4 * math.factorial(l+s) * math.factorial(l-s))
	with sm.evaluate(False):
		Symbolic_sqrt = sm.sqrt( Numerical_part / sm.pi) 

	Suma = 0 
	for r in range(0, l-s +1):
		Suma = Binomial(l-s,r) * Binomial(l+s,r+s-m) * \
			(-1)**( l-r-s) * (sm.cos(theta/2)/sm.sin(theta/2))**(2*r + s - m) \
			+ Suma 

	Suma = sm.simplify((sm.sin(theta/2))**(2*l) * sm.exp(sm.I*m*phi) * Suma)

	power_minus_one= sm.Rational((-1)**(m))
	with sm.evaluate(False):	
		SWSH = power_minus_one * Symbolic_sqrt * Suma 

	return SWSH 

class Yslm:

	def __init__(self, THETA, PHI, s, l, m, lmax, smax):			

		if (lmax < l):
			print("Error!: lmax<l.")
			sys.exit()
		elif(smax < abs(s)):
			print("Error!: smax < s.")
			sys.exit()
		else:
			self.l = l
			self.m = m
			self.s = s
			self.THETA = THETA
			self.PHI   = PHI
			if l < abs(s):
				self.expression = 0
			else:
				self.expression = Formula_of_Yslm(s, l, m)
		
		
	def print_symbolic(self):						

		sm.pprint(self.expression)

	def eval_in_mesh(self):	

		if (self.s == 0 and self.l == 0 and self.m == 0):		
			return be.sqrt(1.0/(4.0*be.pi))*be.ones_like(self.THETA)

		elif(self.expression == 0):

			return 0.*be.ones_like(self.THETA)

		else: 

			theta, phi = sm.symbols('theta phi')
			EXPRE = sm.simplify(self.expression)
			N_SWSH = sm.lambdify( [theta,phi], EXPRE, "numpy")	
			return N_SWSH(self.THETA, self.PHI)